#! /usr/bin/python2.7

## VERSION NOTES
# old version 6, from /Data05/evlee/pipelines/CSAR_structure_hotspots/CSAR_wrapper_bambased_v6.py. Can remove this line in subsequent versions
# new version 1.0.1, updates the locations of called scripts (to make consistent with new structure pipeline) and uses distutils spawn to find non-unix standard but server-wide scripts
# version 1.0.2, now runs + and - strand steps in parallel to save time

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
parser.add_argument('--tag','-t', action='store', help='prefix tag for all output files')
parser.add_argument('--out_dir','-o', help="specify output dir (need path, use trailing '/')")
#features to add later:
#parser.add_argument('--filter_in','-fi', action='store', nargs='?', default='unspecified', help="input file to retain reads in specified regions, should be .bed file of desired regions")
#parser.add_argument('--filter_out','-fo', action='store', nargs='?', default='unspecified', help="input file to filter out reads in specified regions, should be .bed file of unwanted regions")
args=parser.parse_args()

## Check organizational requirements for correct data processing

# # set scripts directory as same directory in which this wrapper is located
# scrDir=os.path.abspath(os.path.join(__file__, os.pardir))

# # check for required scripts
# reqScr=['shuffle_reads_BAM.pl', 'split_coverage_by_chrom.pl', 'generate_empty_coverage_files.rb', 'make_CSAR_files.R', 'run_CSAR_shuffled.R', 'run_CSAR_saturation.R']
# havScr=os.listdir(scrDir)
# getScr=[]
# for needThis in reqScr:
# 	if needThis not in havScr:
# 		getScr.append(needThis)

# if len(getScr) != 0:
# 	for getThis in getScr:
# 		print "(!) Need '%s'" % getThis
# 	raise ValueError("Missing reqired scripts for analysis, move indicated scripts to same directory as this script")

#locations of C, Bash, and R scripts
structure_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
shuffle_reads_BAM=structure_dir+"/CSAR_utils/"+"shuffle_reads_BAM.pl"
split_coverage_by_chrom=structure_dir+"/CSAR_utils/"+"split_coverage_by_chrom.pl"
generate_empty_coverage_files=structure_dir+"/CSAR_utils/"+"generate_empty_coverage_files.rb"
make_CSAR_files=structure_dir+"/CSAR_utils/"+"make_CSAR_files.R"
run_CSAR_shuffled=structure_dir+"/CSAR_utils/"+"run_CSAR_shuffled.R"
run_CSAR_saturation=structure_dir+"/CSAR_utils/"+"run_CSAR_saturation.R"



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

## Begin data processing

# Notify user run is begining
print "Hail CSAR!"


# Shuffle BAM files
print "Shuffling BAM files..."
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rightnow= "_".join(datelist)
dsRNA_BAM_shuffled = args.dsRNA_BAM.replace("bam", rightnow+"_shuffled.bam")
ssRNA_BAM_shuffled = args.ssRNA_BAM.replace("bam", rightnow+"_shuffled.bam")
subprocess.check_call(['perl',shuffle_reads_BAM, args.dsRNA_BAM, args.ssRNA_BAM, dsRNA_BAM_shuffled, ssRNA_BAM_shuffled])


# Calculate genome coverage for shuffled
print "Calculating genome coverage for shuffled BAMs"
for bam in [dsRNA_BAM_shuffled, ssRNA_BAM_shuffled]:
	tag = basename(bam)
	tag = tag.replace("."+rightnow+"_shuffled.bam", "")
	#calculate genome coverage
	plusCov=open(outPrm+'shuffled_'+tag+'.plus.coverage.txt','w')
	genomeCovPlusPs = subprocess.Popen([BEDTOOLS, 'genomecov', '-ibam', bam, '-d', '-split','-strand',  '+'], stdout=plusCov)
	minusCov=open(outPrm+'shuffled_'+tag+'.minus.coverage.txt','w')
	genomeCovMinusPs = subprocess.Popen([BEDTOOLS, 'genomecov', '-ibam', bam, '-d', '-split','-strand',  '-'], stdout=minusCov)
	genomeCovPlusPs.wait()
	genomeCovMinusPs.wait()
	plusCov.close()
	minusCov.close()
	#split by chromosome
	splitPlusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outPrm+'shuffled_'+tag+'.plus.coverage.txt', outPrm+'shuffled_'+tag+'.plus'])
	splitMinusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outPrm+'shuffled_'+tag+'.minus.coverage.txt', outPrm+'shuffled_'+tag+'.minus'])
	splitPlusPs.wait()
	splitMinusPs.wait()
	#remove initial genome coverage files
	subprocess.check_call(['rm', outPrm+'shuffled_'+tag+'.plus.coverage.txt'])
	subprocess.check_call(['rm', outPrm+'shuffled_'+tag+'.minus.coverage.txt'])
	#generate empty coverage files for chromosomes with 0 coverage
	emptyPlusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outPrm+'shuffled_'+tag+'.plus', args.chr_len])
	emptyMinusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outPrm+'shuffled_'+tag+'.minus', args.chr_len])
	emptyPlusPs.wait()
	emptyMinusPs.wait()
	#make CSAR input files
	inputPlusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outPrm, 'shuffled_'+tag+'.plus', 'Forward'])
	inputMinusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outPrm, 'shuffled_'+tag+'.minus', 'Reverse'])
	inputPlusPs.wait()
	inputMinusPs.wait()
	#remove split coverage files
	doneFiles=glob.glob(outPrm+'shuffled_'+tag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outPrm+'shuffled_'+tag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR shuffled
dsRNA_tag = basename(dsRNA_BAM_shuffled)
dsRNA_tag = dsRNA_tag.replace("."+rightnow+"_shuffled.bam", "")
ssRNA_tag = basename(ssRNA_BAM_shuffled)
ssRNA_tag = ssRNA_tag.replace("."+rightnow+"_shuffled.bam", "")
print "Running CSAR for shuffled BAMs..."
subprocess.check_call([RSCRIPT, '--vanilla', run_CSAR_shuffled, outPrm+'shuffled_'+dsRNA_tag, outPrm+'shuffled_'+ssRNA_tag, outPrm+'shuffled_'+args.tag, args.chr_len, outPrm+'shuffled_'+args.tag+'.CSAR.bed', outPrm+'shuffled_'+args.tag+'.CSAR.threshold.txt'])
doneFiles=glob.glob(outPrm+'*CSARNhits')
for target in doneFiles:
	os.remove(target)
doneFiles=glob.glob(outPrm+'*CSARScore')
for target in doneFiles:
	os.remove(target)

# Calculate genome coverage for true
print "Calculating genome coverage for real BAMs"
for bam in [args.dsRNA_BAM, args.ssRNA_BAM]:
	tag = basename(bam)
	tag = tag.replace(".bam", "")
	#calculate genome coverage
	plusCov=open(outDir+tag+'.plus.coverage.txt','w')
	genomeCovPlusPs = subprocess.Popen([BEDTOOLS, 'genomecov', '-ibam', bam, '-d', '-split','-strand',  '+'], stdout=plusCov)
	minusCov=open(outDir+tag+'.minus.coverage.txt','w')
	genomeCovMinusPs = subprocess.Popen([BEDTOOLS, 'genomecov', '-ibam', bam, '-d', '-split','-strand',  '-'], stdout=minusCov)
	genomeCovPlusPs.wait()
	genomeCovMinusPs.wait()
	plusCov.close()
	minusCov.close()
	#split by chromosome
	splitPlusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outDir+tag+'.plus.coverage.txt', outDir+tag+'.plus'])
	splitMinusPs = subprocess.Popen(['perl',split_coverage_by_chrom, outDir+tag+'.minus.coverage.txt', outDir+tag+'.minus'])
	splitPlusPs.wait()
	splitMinusPs.wait()
	#remove initial genome coverage files
	subprocess.check_call(['rm', outDir+tag+'.plus.coverage.txt'])
	subprocess.check_call(['rm', outDir+tag+'.minus.coverage.txt'])
	#generate empty coverage files for chromosomes with 0 coverage
	emptyPlusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outDir+tag+'.plus', args.chr_len])
	emptyMinusPs = subprocess.Popen(['ruby', generate_empty_coverage_files, outDir+tag+'.minus', args.chr_len])
	emptyPlusPs.wait()
	emptyMinusPs.wait()
	#make CSAR input files
	inputPlusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outDir, tag+'.plus', 'Forward'])
	inputMinusPs = subprocess.Popen([RSCRIPT, '--vanilla', make_CSAR_files, outDir, tag+'.minus', 'Reverse'])
	inputPlusPs.wait()
	inputMinusPs.wait()
	#remove split coverage files
	doneFiles=glob.glob(outDir+tag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	doneFiles=glob.glob(outDir+tag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR on true data
dsRNA_tag = basename(args.dsRNA_BAM)
dsRNA_tag = dsRNA_tag.replace(".bam", "")
ssRNA_tag = basename(args.ssRNA_BAM)
ssRNA_tag = ssRNA_tag.replace(".bam", "")
short_ds_tag = dsRNA_tag.split(".")[0]
short_ss_tag = ssRNA_tag.split(".")[0]
subprocess.check_call([RSCRIPT, '--vanilla', run_CSAR_saturation, outDir+dsRNA_tag, outDir+ssRNA_tag, outDir+args.tag, args.chr_len, outPrm+'shuffled_'+args.tag+'.CSAR.threshold.txt', outDir+short_ds_tag+"_enrichedOver_"+short_ss_tag+'.CSAR_peaks.bed', outDir+short_ds_tag+"_enrichedOver_"+short_ss_tag+'.CSAR_peaks_counts.txt'])
doneFiles=glob.glob(outDir+'*CSARNhits')
for target in doneFiles:
	os.remove(target)
doneFiles=glob.glob(outDir+'*CSARScore')
for target in doneFiles:
	os.remove(target)


print "Finished run, output can be found in "+outDir
