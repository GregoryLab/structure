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

#version 1.0


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


#locations of other subprocess scripts that arent present in default PATH
structure_dir=os.path.dirname(os.path.realpath(sys.argv[0]))
BamCoverage=structure_dir+"/"+"BamCoverage"

#command line arguments
parser = argparse.ArgumentParser(description='Calculate basewise structure score over BED intervals using BIGWIG format as an intermediate')
parser.add_argument('in_ref', help='BED12 reference file')
parser.add_argument('in_chr', help='Chromosome length reference file') #Consider revising to extract from BAM header
parser.add_argument ('--multi', help='How to handle coverage of multimappers (w for weight(default), i for ignore multimappings, a for use all mappings, or r for use a random mapping)', default='w')
parser.add_argument('output_folder',help='name of folder to put structure output')
parser.add_argument('output_tag', help='output prefix tag, for files written within output folder')
parser.add_argument('--bam', '-b', action='store', nargs=1, help="mapped reads in BAM format")
args = parser.parse_args()

bam = str(args.bam[0])


#Check for output directory and make it if neccessary
output_folder = re.sub('\/$', '', args.output_folder)
if  os.path.isdir(args.output_folder):
    print "Existing output folder " + output_folder + " detected, will overwrite all internal files"
subprocess.check_call(['mkdir', '-p', output_folder])


# make tmp directory if necessary
tmpDIR=output_folder + '/strand_temp/'
subprocess.check_call(['mkdir', '-p', tmpDIR])

# set "random" tag (unqiue string in case multiple structure runs output to same temp folder)
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rTag= "_".join(datelist)


# Index BAM file (if not already done)
bam_index = bam+".bai" 
if not os.path.isfile(bam_index):
	print >> sys.stderr, "indexing BAM file..."
	indexer = subprocess.Popen(['samtools','index',bam])
	indexer.wait()


# generate bedgraph intermediate files
print >> sys.stderr, "making bedgraphs..."

bgPlus_file = tmpDIR + "tmp_bgP."+rTag+".bgr"
bgMinus_file = tmpDIR + "tmp_bgM."+rTag+".bgr"

script_dir = os.path.dirname(os.path.realpath(__file__))
bgP = subprocess.Popen([BamCoverage ,bam, tmpDIR + "tmp_bgP."+rTag+".bgr","-s+th"+args.multi])
bgM = subprocess.Popen([BamCoverage ,bam, tmpDIR + "tmp_bgP."+rTag+".bgr", "-s-th"+args.multi])

bgP.wait()
bgM.wait()



print >> sys.stderr, "done"


# Generate coverage BIGWIGs
print >> sys.stderr, "making bigwigs..."

bwPlus_file = output_folder + '/' + os.path.basename(bam).replace(".bam", ".plus.bw") 
bwP = subprocess.Popen([BEDGRAPHTOBIGWIG,bgPlus_file,args.in_chr,bwPlus_file])

bwMinus_file = output_folder + '/' + os.path.basename(bam).replace(".bam", ".minus.bw") 
bwM = subprocess.Popen([BEDGRAPHTOBIGWIG,bgMinus_file,args.in_chr,bwMinus_file])

bwP.wait()
bwM.wait()

print >> sys.stderr, "done"


# Paste BIGWIGs into BEDGRAPH
print >> sys.stderr, "pasting bigwigs"

SC_bgCov_file = tmpDIR + "tmp_SC_cov."+rTag+".bgr"
SC_paste = subprocess.Popen([BWTOOL,"paste",bwPlus_file,bwMinus_file,"-o="+SC_bgCov_file])
SC_paste.wait()

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
			ds_cov = set_to_num(DS_old_val) #* DS_norm_factor #remove for strand scores
			ss_cov = set_to_num(SS_old_val) #* SS_norm_factor #remove for strand scores
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

SC_bg_file = tmpDIR + "tmp_SC."+rTag+".bgr"

print >> sys.stderr, "Calculating strand score"

t = multiprocessing.Pool(1)

SC_paste_mod = t.apply_async(struct_thread, args=(SC_bgCov_file, SC_bg_file))

t.close()
t.join()



print >> sys.stderr, "done"


# sort bedgraph output (required for bedGraphToBigWig)
print >> sys.stderr, "Sorting bedgraph output"

os.environ['LC_COLLATE'] = "C" #set case-sensitive sort

SC_bg_sort_file = tmpDIR + "tmp_SC.sorted."+rTag+".bgr"
SC_bg_sort = open(SC_bg_sort_file,'w')
subprocess.check_call(["sort", "-k1,1", "-k2,2n", SC_bg_file], stdout=SC_bg_sort)
SC_bg_sort.close()


print >> sys.stderr, "done"


# Generate strand score BIGWIGs
print >> sys.stderr, "Convert structScore to BIGWIG"
SC_bw_file = output_folder + "/" + args.output_tag + ".strandScore.bw"
subprocess.check_call([BEDGRAPHTOBIGWIG,SC_bg_sort_file,args.in_chr,SC_bw_file])

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

bed_file = tmpDIR+"tmp_bed."+rTag+".txt"
bed = subprocess.Popen([BWTOOL,"extract","bed",tmpRef_file,SC_bw_file,bed_file])
bed.wait()
print >> sys.stderr, "done"


# Create final output by bed sort
print >> sys.stderr, "sorting output"

final_out = open(output_folder + '/'+ args.output_tag + '.strandScore.bed12', 'w')
subprocess.check_call(["sort","-k1,1","-k2,2n","-i",bed_file],stdout=final_out)
final_out.close()

print >> sys.stderr, "done"


# remove temporary directory
subprocess.check_call(["rm","-r",tmpDIR])
