# Structure

Toolkit for nuclease-mediated structure mapping. Includes methods for numerical structure scores, numerical strand bias scores, and structure peaks and valleys.

## Contents

### 1. structScore.py

Calculates basewise structure score over BED intervals using BIGWIG format as an intermediate. Structure scores refer to the log-odds ratio of paired (double-stranded) versus unpaired (single-stranded) conformations. Higher scores imply higher probability of a paired conformation.

### 2. strandScore.py

Calculates basewise strand scores over BED intervals using BIGWIG format as an intermediate. Strand scores refer to the log-odds ratio of sense versus antisense read coverage. Higher scores imply more read bias toward the sense strand. Strand scores are useful for resolving intermolecular duplexes (should contain both sense and antisense reads) from bona-fide intramolecular secondary structure (sense bias). 

### 3. runCSAR.py

Computes structure peaks (highly paired regions) and structure vallyes (highly unpaired regions). Requires lable-shuffled (permutated) datasets for computing a background distribution and empirical false discovery rate (FDR) thresholds.

### 4. extractBed12FromBW.py

Extracts structure or strand scores over BED interviews from BIGWIG files (i.e. those generated with structScore.py or strandScore.py). Saves time when computing scores over new intervals by not repeating whole pipeline. 

### 5. bamRandomMultimapper.py

For each read that maps to multiple loci, extracts one of those loci at random. Replaces 'bam' suffix with 'oneRandomMultimapper.bam' This is one of several ways to treat multi-mapping reads when computing scores or peaks.

### 6. shuffle_reads_BAM.pl

Generate label-shuffled (permuted) datasets in which dsRNA and ssRNA libraries reads are randomly mixed between libraries. This is required for computing a background distribution and empirical false discovery rate (FDR) thresholds. runCSAR.py calls this script within the pipeline, but when performing multiple iterations of peakcalling is generate shuffled readsets only once with shuffle_reads_BAM.pl, and then use the optional shuffled read flags in runCSAR.py.

## Usage

###1. structScore.py

```
python structScore.py --ds_reads [ds_reads] --ss_reads [ss_reads] [in_ref] [in_chr] [output_folder] [output_tag] [OPTIONS]
```

#### Required arguments
	--ds_reads <ds_reads>
		dsRNA (ssRNase-treated) sample mapped reads
	--ss_reads <ss_reads>
		ssRNA (dsRNase-treated) sample mapped reads
	<in_ref>
		BED12 reference file
	<in_chr>
		Chromosome length reference file
	<output_folder>
		name of folder to put structure output
	<output_tag>
		output prefix tag, for files written within output folder

#### Optional arguments:
	--multi <w>      
		How to handle coverage of multimappers 
		(w for weight(default), i for ignore multimappings, a for use all mappings, or r for use a random mapping)


### 2. strandScore.py


```
python strandScore.py --bam [bam] [in_ref] [in_chr] [output_folder] [output_tag] [OPTIONS]
```


#### Required arguments
	--bam, -b  
		mapped reads in BAM format
	<in_ref>
		BED12 reference file
	<in_chr>
		Chromosome length reference file
	<output_folder>
		name of folder to put structure output
	<output_tag>
		output prefix tag, for files written within output folder

#### Optional arguments:
	--multi <w>      
		How to handle coverage of multimappers 
		(w for weight(default), i for ignore multimappings, a for use all mappings, or r for use a random mapping)


###3. runCSAR.py

```
runCSAR.py --dsRNA_BAM [dsRNA_BAM] --ssRNA_BAM [ssRNA_BAM] [chr_len] [OPTIONS]
```

#### Required arguments
	--dsRNA_BAM <dsRNA_BAM>
		dsRNA (ssRNase-treated) sample mapped reads
	--ssRNA_BAM <ssRNA_BAM>
		ssRNA (dsRNase-treated) sample mapped reads
	<chr_len>
		Tab-separated two-column list of chromosome/contig names and lengths

#### Optional arguments:
	--shuffled_dsRNA_BAM, -sds <SHUFFLED_DSRNA_BAM>      
		shuffled dsRNA-seq (i.e. ssRNase-treated) library BAM generated with shuffle_reads_BAM.pl
	--shuffled_ssRNA_BAM, -sss <SHUFFLED_SSRNA_BAM>      
		shuffled ssRNA-seq (i.e. dsRNase-treated) library BAM generated with shuffle_reads_BAM.pl
	--multi <w>      
		How to handle coverage of multimappers 
		(w for weight(default), i for ignore multimappings, a for use all mappings, or r for use a random mapping)
	--tag, -t <TAG>      
		prefix tag for all output files 
	--out_dir, -o <OUT_DIR>      
		specify output dir (need path, use trailing '/')
	--keep_intermediates, -k
		Retain intermediate files	


### 4. extractBed12FromBW.py

```
extractBed12FromBW.py --plus [PLUS] --minus [MINUS] [in_ref] [output_folder] [output_tag]
```

#### Required arguments
	--plus <plus_BIGWIG>
		plus-strand BIGWIG file
	--ssRNA_BAM <minus_BIGWIG>
		minus-strand BIGWIG file
	<in_ref>
		BED12 reference file
	<output_folder>
		name of folder to put structure output
	<output_tag>
		output prefix tag, for files written within output folder


### 5. bamRandomMultimapper.py

```
python bamRandomMultimapper.py [BAM]
```

#### Required arguments
	<bam>
		mapped reads in BAM format


### 6. shuffle_reads_BAM.pl

```
perl shuffle_reads_BAM.pl [sample_fn] [control_fn] [shuffled_sample_fn] [shuffled_control_fn]
```

#### Required arguments
	<sample_fn>
		Sample mapped reads filename (sample vs control is defined arbitrarily)
	<control_fn>
		Control mapped reads filename (sample vs control is defined arbitrarily)
	<shuffled_sample_fn>
		Shuffled mapped reads output filename
	<shuffled_control_fn>
		Shuffled control reads output filename


## Copyright
	Copyright (c) 2013-2018 University of Pennsylvania

	Permission is hereby granted, free of charge, to any person obtaining a
	copy of this software and associated documentation files (the "Software"),
	to deal in the Software without restriction, including without limitation
	the rights to use, copy, modify, merge, publish, distribute, sublicense,
	and/or sell copies of the Software, and to permit persons to whom the
	Software is furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in
	all copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
	OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
	DEALINGS IN THE SOFTWARE.
