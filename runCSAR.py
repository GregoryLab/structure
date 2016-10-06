#! /usr/bin/python2.7

## VERSION NOTES
# old version 6, from /Data05/evlee/pipelines/CSAR_structure_hotspots/CSAR_wrapper_bambased_v6.py. Can remove this line in subsequent versions
# new version 1.0.1, updates the locations of called scripts (to make consistent with new structure pipeline) and uses distutils spawn to find non-unix standard but server-wide scripts
# version 1.0.2, now runs + and - strand steps in parallel to save time
# version 1.0.3, now runs ds and ss steps in parallel to save time
# version 1.1, now allows specification of existing shuffled BAM files
# version 1.2, now includes option to retain intermediate files


# Note that script makes use of a modified ChIPseqScore function in R that removes the strict requirements for writing CSAR output to the same folder as the CSAR script. The function is coded in the the CSAR R scripts, and a reference copy is also in ChIPseqScore_fixed.R 

## Set up python environment
import sys
import os
import argparse
import subprocess
import re
import glob
import random
import datetime
import time
import os.path
from distutils import spawn

RSCRIPT=spawn.find_executable("Rscript")
if RSCRIPT is None:
   print "***ERROR: Rscript is not found"
   sys.exit("Please instal R / Rscript or make sure it is in the PATH")

BEDTOOLS=spawn.find_executable("bedtools")
if BEDTOOLS is None:
   print  "***ERROR: bedtools is not found"
   sys.exit("Please install bedtools or make sure it is in the PATH")


## Define functions
def basename(filepath):
	basename = ""
	if isinstance(filepath, str):
		string = filepath
	else:
		string = str(filepath)
	paths = filepath.split("/")
	basename = paths[len(paths)-1]
	return basename

## Process arguments
parser = argparse.ArgumentParser(description="Input dsRNA-seq and ssRNA-seq files ")
parser.add_argument('chr_len', help="specify path to chromosome length or contig length file")
parser.add_argument('--dsRNA_BAM','-ds', action='store', help='dsRNA-seq (i.e. ssRNase-treated) library BAM')
parser.add_argument('--ssRNA_BAM','-ss', action='store', help='ssRNA-seq (i.e. dsRNase-treated) library BAM')
parser.add_argument('--shuffled_dsRNA_BAM','-sds', action='store', nargs='?', help='shuffled dsRNA-seq (i.e. ssRNase-treated) library BAM (optional)')
parser.add_argument('--shuffled_ssRNA_BAM','-sss', action='store', nargs='?', help='shuffled ssRNA-seq (i.e. dsRNase-treated) library BAM (optional)')
parser.add_argument ('--multi', help='How to handle coverage of multimappers (w for weight(default), i for ignore multimappings, a for use all mappings, or r for use a random mapping)', default='w')
parser.add_argument('--tag','-t', action='store', help='prefix tag for all output files')
parser.add_argument('--out_dir','-o', help="specify output dir (need path, use trailing '/')")
parser.add_argument('--keep_intermediates','-k',action='store_true',help='Retain intermediate files. Useful for debugging.')
#features to add later:
#parser.add_argument('--filter_in','-fi', action='store', nargs='?', default='unspecified', help="input file to retain reads in specified regions, should be .bed file of desired regions")
#parser.add_argument('--filter_out','-fo', action='store', nargs='?', default='unspecified', help="input file to filter out reads in specified regions, should be .bed file of unwanted regions")
args=parser.parse_args()

#locations of C, Bash, and R scripts
structure_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
shuffle_reads_BAM=structure_dir+"/CSAR_utils/"+"shuffle_reads_BAM.pl"
split_coverage_by_chrom=structure_dir+"/CSAR_utils/"+"split_coverage_by_chrom.pl"
generate_empty_coverage_files=structure_dir+"/CSAR_utils/"+"generate_empty_coverage_files.rb"
make_CSAR_files=structure_dir+"/CSAR_utils/"+"make_CSAR_files.R"
run_CSAR_shuffled=structure_dir+"/CSAR_utils/"+"run_CSAR_shuffled.R"
run_CSAR_saturation=structure_dir+"/CSAR_utils/"+"run_CSAR_saturation.R"
BamCoverage=structure_dir+"/"+"BamCoverage"

# generate output directory
if not args.out_dir:
	if not os.path.isdir("./CSAR_output/perms"):
		subprocess.call(["mkdir", "-p", "./CSAR_output/perms"])
	outDir="./CSAR_output/"
	outPrm="./CSAR_output/perms/"
else:
	if args.out_dir[len(args.out_dir)-1] != '/':
		raise NameError("Need to include trailing '/' specified output arg")
	if not os.path.isdir(args.out_dir):
		subprocess.check_call(["mkdir", args.out_dir])
	if not os.path.isdir(args.out_dir+'perms'):
		subprocess.check_call(["mkdir", args.out_dir+'perms'])
	outDir=args.out_dir
	outPrm=args.out_dir+'perms/'

# define unique temp file string
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rightnow= "_".join(datelist)

## Begin data processing

# Notify user run is begining
print "Hail CSAR!"

# Shuffle BAM files
if not args.shuffled_dsRNA_BAM and not args.shuffled_ssRNA_BAM:
	print "Shuffling BAM files..."
	dsRNA_BAM_shuffled = args.dsRNA_BAM.replace("bam", rightnow+"_shuffled.bam")
	ssRNA_BAM_shuffled = args.ssRNA_BAM.replace("bam", rightnow+"_shuffled.bam")
	subprocess.check_call(['perl',shuffle_reads_BAM, args.dsRNA_BAM, args.ssRNA_BAM, dsRNA_BAM_shuffled, ssRNA_BAM_shuffled])

elif args.shuffled_dsRNA_BAM and args.shuffled_ssRNA_BAM:
	print "Using existing shuffled BAM files..."
	dsRNA_BAM_shuffled = args.shuffled_dsRNA_BAM
	ssRNA_BAM_shuffled = args.shuffled_ssRNA_BAM

else:
	raise NameError("Must specify both shuffled_dsRNA_BAM and ssRNA_BAM_shuffled, or specify neither")


####################################################################################################################################################

# Calculate genome coverage for shuffled
print "Calculating genome coverage for shuffled BAMs"
if not args.shuffled_dsRNA_BAM and not args.shuffled_ssRNA_BAM:
	dsTag = basename(dsRNA_BAM_shuffled)
	dsTag = dsTag.replace("."+rightnow+"_shuffled.bam", "")
	ssTag = basename(ssRNA_BAM_shuffled)
	ssTag = ssTag.replace("."+rightnow+"_shuffled.bam", "")
elif args.shuffled_dsRNA_BAM and args.shuffled_ssRNA_BAM:
	dsTagFields = basename(dsRNA_BAM_shuffled).split(".")
	dsTag = ".".join(dsTagFields[:-2])
	ssTagFields = basename(ssRNA_BAM_shuffled).split(".")
	ssTag = ".".join(ssTagFields[:-2])
#calculate genome coverage
dsGenomeCovPlusPs = subprocess.Popen([BamCoverage, dsRNA_BAM_shuffled, outPrm+'shuffled_'+dsTag+'.plus.coverage.txt',"-s+th"+args.multi])
dsGenomeCovMinusPs = subprocess.Popen([BamCoverage, dsRNA_BAM_shuffled, outPrm+'shuffled_'+dsTag+'.minus.coverage.txt',"-s-th"+args.multi])
ssGenomeCovMinusPs = subprocess.Popen([BamCoverage, ssRNA_BAM_shuffled, outPrm+'shuffled_'+ssTag+'.minus.coverage.txt',"-s-th"+args.multi])
ssGenomeCovPlusPs =subprocess.Popen([BamCoverage, ssRNA_BAM_shuffled, outPrm+'shuffled_'+ssTag+'.plus.coverage.txt',"-s+th"+args.multi])

dsGenomeCovPlusPs.wait()
dsGenomeCovMinusPs.wait()
ssGenomeCovPlusPs.wait()
ssGenomeCovMinusPs.wait()




#split by chromosome
dsSplitPlusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outPrm+'shuffled_'+dsTag+'.plus.coverage.txt', outPrm+'shuffled_'+dsTag+'.plus'])
dsSplitMinusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outPrm+'shuffled_'+dsTag+'.minus.coverage.txt', outPrm+'shuffled_'+dsTag+'.minus'])
ssSplitPlusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outPrm+'shuffled_'+ssTag+'.plus.coverage.txt', outPrm+'shuffled_'+ssTag+'.plus'])
ssSplitMinusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outPrm+'shuffled_'+ssTag+'.minus.coverage.txt', outPrm+'shuffled_'+ssTag+'.minus'])
dsSplitPlusPs.wait()
dsSplitMinusPs.wait()
ssSplitPlusPs.wait()
ssSplitMinusPs.wait()
#remove initial genome coverage files
subprocess.check_call(['rm', outPrm+'shuffled_'+dsTag+'.plus.coverage.txt'])
subprocess.check_call(['rm', outPrm+'shuffled_'+dsTag+'.minus.coverage.txt'])
subprocess.check_call(['rm', outPrm+'shuffled_'+ssTag+'.plus.coverage.txt'])
subprocess.check_call(['rm', outPrm+'shuffled_'+ssTag+'.minus.coverage.txt'])
#generate empty coverage files for chromosomes with 0 coverage
dsEmptyPlusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outPrm+'shuffled_'+dsTag+'.plus', args.chr_len])
dsEmptyMinusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outPrm+'shuffled_'+dsTag+'.minus', args.chr_len])
ssEmptyPlusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outPrm+'shuffled_'+ssTag+'.plus', args.chr_len])
ssEmptyMinusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outPrm+'shuffled_'+ssTag+'.minus', args.chr_len])
dsEmptyPlusPs.wait()
dsEmptyMinusPs.wait()
ssEmptyPlusPs.wait()
ssEmptyMinusPs.wait()
#make CSAR input files
dsInputPlusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outPrm, 'shuffled_'+dsTag+'.plus', 'Forward'])
dsInputMinusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outPrm, 'shuffled_'+dsTag+'.minus', 'Reverse'])
ssInputPlusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outPrm, 'shuffled_'+ssTag+'.plus', 'Forward'])
ssInputMinusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outPrm, 'shuffled_'+ssTag+'.minus', 'Reverse'])
dsInputPlusPs.wait()
dsInputMinusPs.wait()
ssInputPlusPs.wait()
ssInputMinusPs.wait()
#remove split coverage files
if not args.keep_intermediates:
	doneFiles=glob.glob(outPrm+'shuffled_'+dsTag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'shuffled_'+dsTag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'shuffled_'+ssTag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'shuffled_'+ssTag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR shuffled
print "Running CSAR for shuffled BAMs..."
subprocess.check_call([RSCRIPT, '--vanilla', run_CSAR_shuffled, outPrm+'shuffled_'+dsTag, outPrm+'shuffled_'+ssTag, outPrm+'shuffled_'+args.tag, args.chr_len, outPrm+'shuffled_'+args.tag+'.CSAR.bed', outPrm+'shuffled_'+args.tag+'.CSAR.threshold.txt'])
if not args.keep_intermediates:
	doneFiles=glob.glob(outPrm+'*CSARNhits')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'*CSARScore')
	for target in doneFiles:
		os.remove(target)

####################################################################################################################################################

## Calculate genome coverage for true
print "Calculating genome coverage for real BAMs"

dsTag = basename(args.dsRNA_BAM)
dsTag = dsTag.replace(".bam", "")
ssTag = basename(args.ssRNA_BAM)
ssTag = ssTag.replace(".bam", "")
#calculate genome coverage
dsGenomeCovPlusPs = subprocess.Popen([BamCoverage, args.dsRNA_BAM, outDir+dsTag+'.plus.coverage.txt',"-s+th"+args.multi])
dsGenomeCovMinusPs = subprocess.Popen([BamCoverage, args.dsRNA_BAM, outDir+dsTag+'.minus.coverage.txt',"-s-th"+args.multi])
ssGenomeCovMinusPs = subprocess.Popen([BamCoverage, args.ssRNA_BAM, outDir+ssTag+'.minus.coverage.txt',"-s-th"+args.multi])
ssGenomeCovPlusPs =subprocess.Popen([BamCoverage, args.ssRNA_BAM, outDir+ssTag+'.plus.coverage.txt',"-s+th"+args.multi])



dsGenomeCovPlusPs.wait()
dsGenomeCovMinusPs.wait()
ssGenomeCovPlusPs.wait()
ssGenomeCovMinusPs.wait()

#split by chromosome
dsSplitPlusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outDir+dsTag+'.plus.coverage.txt', outDir+dsTag+'.plus'])
dsSplitMinusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outDir+dsTag+'.minus.coverage.txt', outDir+dsTag+'.minus'])
ssSplitPlusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outDir+ssTag+'.plus.coverage.txt', outDir+ssTag+'.plus'])
ssSplitMinusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outDir+ssTag+'.minus.coverage.txt', outDir+ssTag+'.minus'])
dsSplitPlusPs.wait()
dsSplitMinusPs.wait()
ssSplitPlusPs.wait()
ssSplitMinusPs.wait()
#remove initial genome coverage files
subprocess.check_call(['rm', outDir+dsTag+'.plus.coverage.txt'])
subprocess.check_call(['rm', outDir+dsTag+'.minus.coverage.txt'])
subprocess.check_call(['rm', outDir+ssTag+'.plus.coverage.txt'])
subprocess.check_call(['rm', outDir+ssTag+'.minus.coverage.txt'])
#generate empty coverage files for chromosomes with 0 coverage
dsEmptyPlusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outDir+dsTag+'.plus', args.chr_len])
dsEmptyMinusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outDir+dsTag+'.minus', args.chr_len])
ssEmptyPlusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outDir+ssTag+'.plus', args.chr_len])
ssEmptyMinusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outDir+ssTag+'.minus', args.chr_len])
dsEmptyPlusPs.wait()
dsEmptyMinusPs.wait()
ssEmptyPlusPs.wait()
ssEmptyMinusPs.wait()
#make CSAR input files
dsInputPlusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outDir, dsTag+'.plus', 'Forward'])
dsInputMinusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outDir, dsTag+'.minus', 'Reverse'])
ssInputPlusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outDir, ssTag+'.plus', 'Forward'])
ssInputMinusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outDir, ssTag+'.minus', 'Reverse'])
dsInputPlusPs.wait()
dsInputMinusPs.wait()
ssInputPlusPs.wait()
ssInputMinusPs.wait()
#remove split coverage files
if not args.keep_intermediates:
	doneFiles=glob.glob(outDir+dsTag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+dsTag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+ssTag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+ssTag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR on true data
dsTag = basename(args.dsRNA_BAM)
dsTag = dsTag.replace(".bam", "")
ssTag = basename(args.ssRNA_BAM)
ssTag = ssTag.replace(".bam", "")
short_ds_tag = dsTag.split(".")[0]
short_ss_tag = ssTag.split(".")[0]
subprocess.check_call([RSCRIPT, '--vanilla', run_CSAR_saturation, outDir+dsTag, outDir+ssTag, outDir+args.tag, args.chr_len, outPrm+'shuffled_'+args.tag+'.CSAR.threshold.txt', outDir+short_ds_tag+"_enrichedOver_"+short_ss_tag+'.CSAR_peaks.bed', outDir+short_ds_tag+"_enrichedOver_"+short_ss_tag+'.CSAR_peaks_counts.txt'])
if not args.keep_intermediates:
	doneFiles=glob.glob(outDir+'*CSARNhits')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+'*CSARScore')
	for target in doneFiles:
		os.remove(target)

####################################################################################################################################################

print "Finished run, output can be found in "+outDir
