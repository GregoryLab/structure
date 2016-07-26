#! /usr/bin/python2.7
# Wrapper will take .bam alignments for a PIP-seq tetrad
# must have .bam as follows ./[project_dir]/data/topOut/tophat_output_RD_[sample_id]/NR_hits.bam
# must have scripts in ./[project_dir]/scripts/
# run wrapper in ./[project_dir]/data/topOut
# will generate directories (default: CSAR_output/perms) in ./[project_dir]/data/
# use a tag shared by only the 4 PIP-seq samples to run
# Follow file naming format:'[unique_tag]_(d/s)sRNase_(no/)protein_((un/)trimmed/merged)/NR_hits.bam'
# Examples: 'Pp_dark_ssRNase_protein_trimmed/NR_hits.bam', 'totalArab_nuc_rep4_dsRNase_noprotein_merged/NR_hits.bam'
# Ideally, you will use this in conjunction with 'alignment_wrapper.py', which will reduce naming issues
# This is a wrapper, so you need the other scripts to actually do the work. Here's the list:
# bedtools (suite), 'shuffle_reads_BAM.pl', 'split_coverage_by_chrom.pl', 'generate_empty_coverage_files.rb', 'make_CSAR_files.R', 'run_CSAR_shuffled.R', 'run_CSAR_saturation.R'
# all scripts must be in same dir as wrapper script
# !!! DO NOT perform runs with identical 'alignment_tag' but different 'trimming_status'
# This wrapper is an organized and slightly edited version of the PIP-seq analysis piple developed by Fan Li

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
#parser.add_argument('--filter_in','-fi', action='store', nargs='?', default='unspecified', help="input file to retain reads in specified regions, should be .bed file of desired regions")
#parser.add_argument('--filter_out','-fo', action='store', nargs='?', default='unspecified', help="input file to filter out reads in specified regions, should be .bed file of unwanted regions")
args=parser.parse_args()

## Check organizational requirements for correct data processing

# set scripts directory as same directory in which this wrapper is located
scrDir=os.path.abspath(os.path.join(__file__, os.pardir))

# check for required scripts
reqScr=['shuffle_reads_BAM.pl', 'split_coverage_by_chrom.pl', 'generate_empty_coverage_files.rb', 'make_CSAR_files.R', 'run_CSAR_shuffled.R', 'run_CSAR_saturation.R']
havScr=os.listdir(scrDir)
getScr=[]
for needThis in reqScr:
	if needThis not in havScr:
		getScr.append(needThis)

if len(getScr) != 0:
	for getThis in getScr:
		print "(!) Need '%s'" % getThis
	raise ValueError("Missing reqired scripts for analysis, move indicated scripts to same directory as this script")

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


# # set target .bam files and filter if neccissary
# if args.filter_bam:
# 		rTag=str(random.randint(10000,99999))
# 		filtered_bam=open(topDir+"/tmp"+rTag+".bam",'w')
# 		subprocess.check_call(['intersectBed','-abam',topDir+'/NR_hits.bam','-b',args.filter_bam,'-s','-v'],stdout=filtered_bam)
# 		filtered_bam.close()
# 		subprocess.check_call(['samtools','sort',topDir+"/tmp"+rTag+".bam",topDir+'/'+set_bam])
# 		os.remove(topDir+"/tmp"+rTag+".bam")

# Shuffle BAM files
print "Shuffling BAM files..."
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rightnow= "_".join(datelist)
dsRNA_BAM_shuffled = args.dsRNA_BAM.replace("bam", rightnow+"_shuffled.bam")
ssRNA_BAM_shuffled = args.ssRNA_BAM.replace("bam", rightnow+"_shuffled.bam")
#subprocess.check_call(['perl',scrDir+"/"+'shuffle_reads_BAM.pl', args.dsRNA_BAM, args.ssRNA_BAM, dsRNA_BAM_shuffled, ssRNA_BAM_shuffled])


# Calculate genome coverage for shuffled
print "Calculating genome coverage for shuffled BAMs"
for bam in [dsRNA_BAM_shuffled, ssRNA_BAM_shuffled]:
	tag = basename(bam)
	tag = tag.replace(".shuffled.bam", "")
	plusCov=open(outPrm+'shuffled_'+tag+'.plus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', bam, '-d', '-split','-strand',  '+'], stdout=plusCov)
	plusCov.close()
	subprocess.check_call(['perl',scrDir+"/"+'split_coverage_by_chrom.pl', outPrm+'shuffled_'+tag+'.plus.coverage.txt', outPrm+'shuffled_'+tag+'.plus'])
	subprocess.check_call(['rm', outPrm+'shuffled_'+tag+'.plus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+"/"+'generate_empty_coverage_files.rb', outPrm+'shuffled_'+tag+'.plus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+"/"+'make_CSAR_files.R', outPrm, 'shuffled_'+tag+'.plus', 'Forward'])
	doneFiles=glob.glob(outPrm+'shuffled_'+tag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	minusCov=open(outPrm+'shuffled_'+tag+'.minus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', bam, '-d', '-split','-strand',  '-'], stdout=minusCov)
	minusCov.close()
	subprocess.check_call(['perl',scrDir+"/"+'split_coverage_by_chrom.pl', outPrm+'shuffled_'+tag+'.minus.coverage.txt', outPrm+'shuffled_'+tag+'.minus'])
	subprocess.check_call(['rm', outPrm+'shuffled_'+tag+'.minus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+"/"+'generate_empty_coverage_files.rb', outPrm+'shuffled_'+tag+'.minus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+"/"+'make_CSAR_files.R', outPrm, 'shuffled_'+tag+'.minus', 'Reverse'])
	doneFiles=glob.glob(outPrm+'shuffled_'+tag+'.minus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)

# Run CSAR shuffled
dsRNA_tag = basename(dsRNA_BAM_shuffled)
dsRNA_tag = dsRNA_tag.replace("."+rightnow+"_shuffled.bam", "")
ssRNA_tag = basename(ssRNA_BAM_shuffled)
ssRNA_tag = ssRNA_tag.replace("."+rightnow+"_shuffled.bam", "")
print "Running CSAR for shuffled BAMs..."
subprocess.check_call(['/usr/bin/Rscript', '--vanilla', scrDir+"/"+'run_CSAR_shuffled.R', outPrm+'shuffled_'+dsRNA_tag, outPrm+'shuffled_'+ssRNA_tag, outPrm+'shuffled_'+args.tag, args.chr_len, outPrm+'shuffled_'+args.tag+'.CSAR.bed', outPrm+'shuffled_'+args.tag+'.CSAR.threshold.txt'])
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
	plusCov=open(outDir+tag+'.plus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', bam, '-d', '-split', '-strand', '+'], stdout=plusCov)
	plusCov.close()
	subprocess.check_call(['perl',scrDir+"/"+'split_coverage_by_chrom.pl', outDir+tag+'.plus.coverage.txt', outDir+tag+'.plus'])
	subprocess.check_call(['rm', outDir+tag+'.plus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+"/"+'generate_empty_coverage_files.rb', outDir+tag+'.plus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+"/"+'make_CSAR_files.R', outDir, tag+'.plus', 'Forward'])
	doneFiles=glob.glob(outDir+tag+'.plus.*.coverage.txt')
	for target in doneFiles:
		os.remove(target)
	minusCov=open(outDir+tag+'.minus.coverage.txt','w')
	subprocess.check_call(['genomeCoverageBed', '-ibam', bam, '-d', '-split', '-strand', '-'], stdout=minusCov)
	minusCov.close()
	subprocess.check_call(['perl',scrDir+"/"+'split_coverage_by_chrom.pl', outDir+tag+'.minus.coverage.txt', outDir+tag+'.minus'])
	subprocess.check_call(['rm', outDir+tag+'.minus.coverage.txt'])
	subprocess.check_call(['ruby', scrDir+"/"+'generate_empty_coverage_files.rb', outDir+tag+'.minus', args.chr_len])
	subprocess.check_call(['Rscript', '--vanilla', scrDir+"/"+'make_CSAR_files.R', outDir, tag+'.minus', 'Reverse'])
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
subprocess.check_call(['/usr/bin/Rscript', '--vanilla', scrDir+"/"+'run_CSAR_saturation.R', outDir+dsRNA_tag, outDir+ssRNA_tag, outDir+args.tag, args.chr_len, outPrm+'shuffled_'+args.tag+'.CSAR.threshold.txt', outDir+short_ds_tag+"_enrichedOver_"+short_ss_tag+'.CSAR_peaks.bed', outDir+short_ds_tag+"_enrichedOver_"+short_ss_tag+'.CSAR_peaks_counts.txt'])
doneFiles=glob.glob(outDir+'*CSARNhits')
for target in doneFiles:
	os.remove(target)
doneFiles=glob.glob(outDir+'*CSARScore')
for target in doneFiles:
	os.remove(target)


print "Finished run, output can be found in "+outDir