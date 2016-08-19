#!/usr/local/bin/python
#
#  Copyright (c) 2016 University of Pennsylvania
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.

#old version 7 from /Data05/evlee/utils folder (i.e. scripts based on those from Sager). Can remove this line in subsequent versions
#new version 1.0


import sys
import os
import argparse
import subprocess
import re
import random
import math
import multiprocessing
import datetime
from distutils import spawn


#check python version
if sys.hexversion < 0x020700F0:
    print "Detected Python " + sys.version
    sys.exit("***ERROR: Must be using Python 2.7.x (recommended) or above")


#find subprocess scripts that aren't standard in unix environement
BWTOOL=spawn.find_executable("bwtool")
if BWTOOL is None:
   print "***ERROR: bwtool is not found"
   sys.exit("Please instal bwtool or make sure it is in the PATH")

BEDTOOLS=spawn.find_executable("bedtools")
if BEDTOOLS is None:
   print  "***ERROR: bedtools is not found"
   sys.exit("Please install bedtools or make sure it is in the PATH")

SAMTOOLS=spawn.find_executable("samtools")
if SAMTOOLS is None:
   print  "***ERROR: samtools is not found"
   sys.exit("Please install samtools or make sure it is in the PATH")

BEDGRAPHTOBIGWIG=spawn.find_executable("bedGraphToBigWig")
if BEDGRAPHTOBIGWIG is None:
   print  "***ERROR: bedGraphToBigWig is not found"
   sys.exit("Please install bedGraphToBigWig or make sure it is in the PATH")


#command line arguments
parser = argparse.ArgumentParser(description='Calculate basewise structure score over BED intervals using BIGWIG format as an intermediate')
parser.add_argument('in_ref', help='BED12 reference file')
parser.add_argument('in_chr', help='Chromosome length reference file') #Consider revising to use standard format, e.g. .dict files from GATK toolkit
parser.add_argument('output_folder',help='name of folder to put structure output')
parser.add_argument('output_tag', help='BED6 files with comma seperated list of scores in columns 7+')
parser.add_argument ('--multi', help='How to handle coverage of multimappers (w for weight(default), i for ignore multimappings, a for use all mappings, or r for use a random mapping)', default='w')
parser.add_argument('--ds_reads', '-ds', action='store', nargs=1, help="dsRNA-seq reads in BAM format, this library was generated using ssRNase (RNaseONE) to isolate dsRNA fragments")
parser.add_argument('--ss_reads', '-ss', action='store', nargs=1, help="ssRNA-seq reads in BAM format, this library was generated using dsRNase (RNaseV1) to isolate ssRNA fragments")
args = parser.parse_args()

ds_reads = str(args.ds_reads[0])
ss_reads = str(args.ss_reads[0])


#Check for output directory and make it if neccessary
output_folder = re.sub('\/$', '', args.output_folder)
if  os.path.isdir(args.output_folder):
    print "Existing output folder " + output_folder + " detected, will overwrite all internal files"
subprocess.check_call(['mkdir', '-p', output_folder])


# make tmp directory if necessary
tmpDIR=output_folder + '/structure_temp/'
subprocess.check_call(['mkdir', '-p', tmpDIR])

# set "random" tag (unqiue string in case multiple structure runs output to same temp folder)
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rTag= "_".join(datelist)


# Index BAM files (if not already done)
ds_reads_index = ds_reads+".bai" 
ss_reads_index = ss_reads+".bai" 
if not os.path.isfile(ds_reads_index):
	print >> sys.stderr, "indexing dsRNA-seq BAM file..."
	DS_indexer = subprocess.Popen(['samtools','index',ds_reads])
	DS_indexer.wait()
if not os.path.isfile(ss_reads_index):
	print >> sys.stderr, "indexing ssRNA-seq BAM file..."
	SS_indexer = subprocess.Popen(['samtools','index',ss_reads])
	SS_indexer.wait()


#Define normalization factors by counting number of matched BASES in primary read alignments (i.e. number of unique read IDs)
#0x100 flag is secondary alignment in SAM v1 specifications, -F filters all reads with this flag
print >> sys.stderr, "counting number of matched positions in primary read alignments..."
match_grab = ur"([0-9]+)M"

DS_mapped_reads = 0
DS_read_dist = 0
ds_uniq_readID = subprocess.Popen([SAMTOOLS, 'view', '-F0x100',ds_reads], stdout=subprocess.PIPE)
for line in ds_uniq_readID.stdout:
	fields = line.rstrip().split("\t")
	CIGAR = fields[5] #based on SAM v1 specifications
	DS_mapped_reads += 1
	match_len = re.findall(match_grab,CIGAR)
	match_len = [ int(x) for x in match_len ]
	DS_read_dist += sum(match_len)
ds_uniq_readID.stdout.close()	
	
SS_mapped_reads = 0
SS_read_dist = 0
ss_uniq_readID = subprocess.Popen([SAMTOOLS, 'view', '-F0x100',ss_reads], stdout=subprocess.PIPE)
for line in ss_uniq_readID.stdout:
	fields = line.rstrip().split("\t")
	CIGAR = fields[5] #based on SAM v1 specifications
	SS_mapped_reads += 1
	match_len = re.findall(match_grab,CIGAR)
	match_len = [ int(x) for x in match_len ]
	SS_read_dist += sum(match_len)
ss_uniq_readID.stdout.close()	

print >> sys.stderr, "calculating normalization factors..."

DS_avg_read_len = float(DS_read_dist) / DS_mapped_reads
SS_avg_read_len = float(SS_read_dist) / SS_mapped_reads
base_factor = float(max(DS_read_dist,SS_read_dist))
DS_norm_factor = base_factor / DS_read_dist
SS_norm_factor = base_factor / SS_read_dist

print >> sys.stderr, "Unique mapped read IDs, dsRNA-seq: " + str(DS_mapped_reads)
print >> sys.stderr, "Avg mapped read length, dsRNA-seq: " + str(DS_avg_read_len)
print >> sys.stderr, "Total matched read len, dsRNA-seq: " + str(DS_read_dist)
print >> sys.stderr, "Read normalization factor, dsRNA-seq: " + str(DS_norm_factor)

print >> sys.stderr, "Unique mapped read IDs, ssRNA-seq: " + str(SS_mapped_reads)
print >> sys.stderr, "Avg mapped read length, ssRNA-seq: " + str(SS_avg_read_len)
print >> sys.stderr, "Total matched read len, ssRNA-seq: " + str(SS_read_dist)
print >> sys.stderr, "Read normalization factor, ssRNA-seq: " + str(SS_norm_factor)


# generate bedgraph intermediate files
print >> sys.stderr, "making bedgraphs..."

DS_bgPlus_file = tmpDIR + "tmp_DS_bgP."+rTag+".bgr"
DS_bgMinus_file = tmpDIR + "tmp_DS_bgM."+rTag+".bgr"
SS_bgPlus_file = tmpDIR + "tmp_SS_bgP."+rTag+".bgr"
SS_bgMinus_file = tmpDIR + "tmp_SS_bgM."+rTag+".bgr"


DS_bgP =  subprocess.Popen(["./BamCoverage",bam,DS_bgPlus_file,"-s+t"+args.multi])
DS_bgM = subprocess.Popen(["./BamCoverage",bam,DS_bgMinus_file,"-s+t"+args.multi])
SS_bgP = subprocess.Popen(["./BamCoverage",bam,SS_bgPlus_file,"-s+t"+args.multi])
SS_bgM = subprocess.Popen(["./BamCoverage",bam,SS_bgMinus_file,"-s+t"+args.multi])


DS_bgP.wait()
DS_bgM.wait()
SS_bgP.wait()
SS_bgM.wait()

print >> sys.stderr, "done"



# Generate coverage BIGWIGs
print >> sys.stderr, "making bigwigs..."

DS_bwPlus_file = output_folder + '/' + os.path.basename(ds_reads).replace(".bam", ".plus.bw") 
DS_bwP = subprocess.Popen([BEDGRAPHTOBIGWIG,DS_bgPlus_file,args.in_chr,DS_bwPlus_file])

DS_bwMinus_file = output_folder + '/' + os.path.basename(ds_reads).replace(".bam", ".minus.bw") 
DS_bwM = subprocess.Popen([BEDGRAPHTOBIGWIG,DS_bgMinus_file,args.in_chr,DS_bwMinus_file])

SS_bwPlus_file = output_folder + '/' + os.path.basename(ss_reads).replace(".bam", ".plus.bw") 
SS_bwP = subprocess.Popen([BEDGRAPHTOBIGWIG,SS_bgPlus_file,args.in_chr,SS_bwPlus_file])

SS_bwMinus_file = output_folder + '/' + os.path.basename(ss_reads).replace(".bam", ".minus.bw")
SS_bwM = subprocess.Popen([BEDGRAPHTOBIGWIG,SS_bgMinus_file,args.in_chr,SS_bwMinus_file])

DS_bwP.wait()
DS_bwM.wait()
SS_bwP.wait()
SS_bwM.wait()
print >> sys.stderr, "done"


# Paste BIGWIGs into BEDGRAPH
print >> sys.stderr, "pasting bigwigs"

SC_bgCovP_file = tmpDIR + "tmp_SC_covP."+rTag+".bgr"
SC_pasteP = subprocess.Popen([BWTOOL,"paste",DS_bwPlus_file,SS_bwPlus_file,"-o="+SC_bgCovP_file])
SC_bgCovM_file = tmpDIR + "tmp_SC_covM."+rTag+".bgr"
SC_pasteM = subprocess.Popen([BWTOOL,"paste",DS_bwMinus_file,SS_bwMinus_file,"-o="+SC_bgCovM_file])

SC_pasteP.wait()
SC_pasteM.wait()

print >> sys.stderr, "done"


# Calculate structure scores genome wide
def glog2(x):
	y = x + math.sqrt( 1 + (x*x) )
	return math.log(y,2)

def set_to_num(x):
	try:
		return float(x)
	except ValueError:
		return 0
def calc_structScore(input_line,start_set,stop_set,old_chr,DS_old_val,SS_old_val):
	out_line = "hold"
	(chrom,start,stop,DS_val,SS_val) = input_line.rstrip().split("\t")
	if DS_val == DS_old_val and SS_val == SS_old_val and chrom == old_chr:
		stop_set = stop
	else:
		if old_chr != '' and ( DS_old_val != "NA" or SS_old_val != "NA"):
			ds_cov = set_to_num(DS_old_val) * DS_norm_factor
			ss_cov = set_to_num(SS_old_val) * SS_norm_factor
			structScore = glog2(ds_cov) - glog2(ss_cov)
			out_line = "\t".join([old_chr,start_set,stop_set,str(structScore)]) + "\n"
		start_set = start
		stop_set = stop
	old_chr = chrom
	DS_old_val = DS_val
	SS_old_val = SS_val
	return [out_line,start_set,stop_set,old_chr,DS_old_val,SS_old_val]

def struct_thread(input_file, output_file):
	try:
		input_pipe = open(input_file, 'r')
		output_pipe = open(output_file, 'w')
		DS_old_val = ''
		SS_old_val = ''
		old_chr = ''
		start_set = ''
		stop_set = ''
#		for line in iter(input_pipe.readline, ''):
		for line in input_pipe:
			(chrom,start,stop,DS_val,SS_val) = line.rstrip().split("\t")
			(out_line,start_set,stop_set,old_chr,DS_old_val,SS_old_val) = calc_structScore(line,start_set,stop_set,old_chr,DS_old_val,SS_old_val)
			if out_line != "hold":
				output_pipe.write(out_line)
	finally:
		try:
			output_pipe.close()
		finally:
			input_pipe.close()

SC_bgP_file = tmpDIR + "tmp_SC_P."+rTag+".bgr"
SC_bgM_file = tmpDIR + "tmp_SC_M."+rTag+".bgr"

print >> sys.stderr, "Calculating structure score"

t = multiprocessing.Pool(2)

SC_pasteP_mod = t.apply_async(struct_thread, args=(SC_bgCovP_file, SC_bgP_file))
SC_pasteM_mod = t.apply_async(struct_thread, args=(SC_bgCovM_file, SC_bgM_file))

t.close()
t.join()

print >> sys.stderr, "done"


# sort bedgraph output (required for bedGraphToBigWig)
print >> sys.stderr, "Sorting bedgraph output"

os.environ['LC_COLLATE'] = "C" #set case-sensitive sort

SC_bgP_sort_file = tmpDIR + "tmp_SC_P.sorted."+rTag+".bgr"
SC_bgP_sort = open(SC_bgP_sort_file,'w')
subprocess.check_call(["sort", "-k1,1", "-k2,2n", SC_bgP_file], stdout=SC_bgP_sort)
SC_bgP_sort.close()

SC_bgM_sort_file = tmpDIR + "tmp_SC_M.sorted."+rTag+".bgr"
SC_bgM_sort = open(SC_bgM_sort_file,'w')
subprocess.check_call(["sort", "-k1,1", "-k2,2n", SC_bgM_file], stdout=SC_bgM_sort)
SC_bgM_sort.close()

print >> sys.stderr, "done"


# Generate structure score BIGWIGs
print >> sys.stderr, "Convert structScore to BIGWIG"
SC_bwP_file = output_folder + "/" + args.output_tag + ".plus.bw"
SC_mkbwP = subprocess.Popen([BEDGRAPHTOBIGWIG,SC_bgP_sort_file,args.in_chr,SC_bwP_file])

SC_bwM_file = output_folder + "/" + args.output_tag + ".minus.bw"
SC_mkbwM = subprocess.Popen([BEDGRAPHTOBIGWIG,SC_bgM_sort_file,args.in_chr,SC_bwM_file])
SC_mkbwP.wait()
SC_mkbwM.wait()
print >> sys.stderr, "done"


# reformat bed file reference
openRef = open(args.in_ref,'r')
tmpRef_file = tmpDIR+"tmp_ref."+rTag+".txt"
tmpRef_open = open(tmpRef_file,'w')

print >> sys.stderr, "reformating reference"
for line in openRef:
	j = line.rstrip().split("\t")
	tmpRef_open.write("\t".join([j[0],j[1],j[2],":".join(j[3:len(j)])]) + "\n")

print >> sys.stderr, "done"

openRef.close()
tmpRef_open.close()


# transfer data to BED reference
print >> sys.stderr, "extracting coverage to BED format"

bedPlus_file = tmpDIR+"tmp_bedP."+rTag+".txt"
bedP = subprocess.Popen([BWTOOL,"extract","bed",tmpRef_file,SC_bwP_file,bedPlus_file])
bedMinus_file = tmpDIR+"tmp_bedM."+rTag+".txt"
bedM = subprocess.Popen([BWTOOL,"extract","bed",tmpRef_file,SC_bwM_file,bedMinus_file])
bedP.wait()
bedM.wait()
print >> sys.stderr, "done"


# filter strand data
bedAll_file = tmpDIR+"tmp_bedA."+rTag+".txt"
bedAll_open = open(bedAll_file, 'w')

bedPlus_open = open(bedPlus_file,'r')
bedMinus_open = open(bedMinus_file,'r')

print >> sys.stderr, "filtering strand data"

for line in bedPlus_open:
	j = line.rstrip().split(":")
	k = "\t".join(j)
	l = k.split("\t") #splitting based on ':' and '\t'?
	go_line = l[:len(l)-2] + l[len(l)-1:]
	if j[2] == "+":
		bedAll_open.write("\t".join(go_line) + "\n")

bedPlus_open.close()

for line in bedMinus_open:
	j = line.rstrip().split(":")
	k = "\t".join(j)
	l = k.split("\t")
	go_line = l[:len(l)-2] + l[len(l)-1:]
	if j[2] == "-":
		bedAll_open.write("\t".join(go_line) + "\n")

bedMinus_open.close()
bedAll_open.close()

print >> sys.stderr, "done"


# Create final output by bed sort
print >> sys.stderr, "sorting output"

final_out = open(output_folder + '/'+ args.output_tag + '.bed12', 'w')
subprocess.check_call(["sort","-k1,1","-k2,2n","-i",bedAll_file],stdout=final_out)
final_out.close()

print >> sys.stderr, "done"


# remove temporary directory
subprocess.check_call(["rm","-r",tmpDIR])
