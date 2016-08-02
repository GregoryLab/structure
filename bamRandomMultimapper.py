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
import pysam
from distutils import spawn


#check python version
if sys.hexversion < 0x020700F0:
    print "Detected Python " + sys.version
    sys.exit("***ERROR: Must be using Python 2.7.x (recommended) or above")

#find subprocess scripts that aren't standard in unix environement
SAMTOOLS=spawn.find_executable("samtools")
if SAMTOOLS is None:
   print  "***ERROR: samtools is not found"
   sys.exit("Please install samtools or make sure it is in the PATH")

#command line arguments
parser = argparse.ArgumentParser(description='Calculate basewise structure score over BED intervals using BIGWIG format as an intermediate')
parser.add_argument('input_bam', help='BAM file with multimappers')
args = parser.parse_args()

# set "random" tag (unqiue string in case multiple structure runs output to same temp folder)
now = datetime.datetime.now()
datelist = [str(now.year),str(now.month),str(now.day),str(now.hour),str(now.minute),str(now.second),str(now.microsecond)]
rTag= "_".join(datelist)

# sort BAM by name
print "sorting input BAM by name..."
bamNameSorted = args.input_bam.replace("bam", "")+rTag+".nameSorted.bam"
subprocess.check_call([SAMTOOLS, "sort", "-n", args.input_bam, "-f", bamNameSorted])

# set random seed (for replicable behavior when rerunning script)
random.seed(18)

## Run through each line of bam file, defining read name groups. For each read name group, pull a random read
print "pulling random instance of each multimapping read..."
bamFile = pysam.AlignmentFile(bamNameSorted, "rb")
outputBamFile = str(args.input_bam).replace("bam", "oneRandomMultimapper.nameSorted.bam")
outputBamFileFH = pysam.AlignmentFile(outputBamFile, "wb", template=bamFile)


current_id = ""
reads = []
for read in bamFile:
	if current_id == "": #initialize with first line
		current_id = str(read.query_name)
		reads.append(read)
	else:
		read_id = str(read.query_name)
		if read_id == current_id:
			reads.append(read)
		else:
			output_read = random.choice(reads)
			outputBamFileFH.write(output_read)
			current_id = str(read_id)
			reads = []
			reads.append(read)

#print out final read group
output_read = random.choice(reads)
outputBamFileFH.write(output_read)

bamFile.close()
outputBamFileFH.close()

# sort output BAM by position
print "sorting output BAM by position..."
outputBamFileSorted = str(args.input_bam).replace("bam", "oneRandomMultimapper.bam")
subprocess.check_call([SAMTOOLS, "sort", outputBamFile, "-f", outputBamFileSorted])

#remove intermediate nameSorted files
os.remove(outputBamFile)
os.remove(bamNameSorted)
