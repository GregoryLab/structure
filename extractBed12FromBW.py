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

###extractBed12FromBW.py###

#Use for reanalyzing structure or strand scores with a new set of bed12-formatted genomic intervals (e.g. transcriptome)

#version1.0, based on the last steps of structScore.py.


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


#find subprocess scripts that aren't standard in unix environment
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
parser.add_argument('output_folder',help='name of folder to put output files')
parser.add_argument('output_tag', help='tag (prefix) for output files')
parser.add_argument('--plus', '-p', action='store', nargs=1, help="plus strand BigWig file")
parser.add_argument('--minus', '-m', action='store', nargs=1, help="minus strand BiwWig file")
args = parser.parse_args()

plus = str(args.minus[0])
minus = str(args.plus[0])


#Check for output directory and make it if neccessary
output_folder = re.sub('\/$', '', args.output_folder)
if  os.path.isdir(args.output_folder):
    print "Existing output folder " + output_folder + " detected, will overwrite internal files with same name as output files"
subprocess.check_call(['mkdir', '-p', output_folder])


# make tmp directory if necessary
tmpDIR=output_folder + '/structure_temp/'
subprocess.check_call(['mkdir', '-p', tmpDIR])

# set "random" tag (unqiue string in case multiple structure runs output to same temp folder)
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rTag= "_".join(datelist)


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
bedP = subprocess.Popen([BWTOOL,"extract","bed",tmpRef_file,plus,bedPlus_file])
bedMinus_file = tmpDIR+"tmp_bedM."+rTag+".txt"
bedM = subprocess.Popen([BWTOOL,"extract","bed",tmpRef_file,minus,bedMinus_file])
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
