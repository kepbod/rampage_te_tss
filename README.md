# Pipeline to identify expressed TE-derived TSSs using RAMPAGE

## A schematic flow shows the pipeline

![workflow](https://github.com/kepbod/rampage_te_tss/blob/master/workflow.jpg)

## Prerequisites

### Softwares

* [StringTie](https://ccb.jhu.edu/software/stringtie/)
* [F-seq](http://fureylab.web.unc.edu/software/fseq/)

### Python libraries 

* [docopt](http://docopt.org/)
* [seqlib](https://github.com/kepbod/seqlib)
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
* [pybedtools](https://daler.github.io/pybedtools/)
* [joblib](https://joblib.readthedocs.io/en/latest/)

## Usage

### Step 1:

Fetch proper read pairs and remove PCR dunplicates.

```
Usage: rm_pcr.py [options] <rampage>...

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output directory. [default: rampage_peak]
    --min=MIN                      Minimum read counts. [default: 1]
    --uniq                         Retain read pairs with unique mate2.
```

* Inputs: BAM files of RAMPAGE (`<rampage>...`)
* Output: A output folder containing relevant files (`-o OUTPUT`)
    * `rampage_plus_5end.bed`: BED file of the 5' end of plus strand read pairs
    * `rampage_plus_3read.bed`: BED file of the 3' end of plus strand read pairs
    * `rampage_minus_5end.bed`: BED file of the 5' end of minus strand read pairs
    * `rampage_minus_3read.bed`: BED file of the 3' end of minus strand read pairs
    * `rampage_link.bed`: BED file linking the 5' and 3' ends of read pairs

Note:
1. If there are multiple RAMPAGE BAM files (different replicates) derived from the same samples, you could simply list them afterwards.
2. You could run with multiple threads using `-p THREAD`.
3. You could set the minimum read filter for read pairs using `--min=MIN`.
4. You could set `--uniq` to retain read pairs with unique mate2.

Example: 

```
rm_pcr.py -o rampage_TETSS_peak --uniq rampage_rep1.bam rampage_rep2.bam
```

### Step 2:

Call peaks using 5' end of RAMPAGE read pairs

```
Usage: call_peak.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
    --wig                          Create Wig files.
    -p PERCENT                     Retained percent of reads in resized peaks.
                                   [default: 0.7]
```

* Input: the output folder created by `rm_pcr.py`
* Output: `rampage_peaks.txt` under the input folder

Format of `rampage_peaks.txt`:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome                    |
| Start       | Start of peak region          |
| End         | End of peak region            |
| Name        | peak                          |
| Score       | 0                             |
| Strand      | Strand of peak                |
| Peak        | peak site                     |
| Height      | Height of peak site           |
| Peak_reads  | Reads of (peak site ± 2 bp)   |
| Total       | Total reads of peak region    |
| Start_Fseq  | Start of F-seq peak region    |
| End_Fseq    | End of F-seq peak region      |
| RPM         | RPM of peak region            |

Note:

1. You could set feature length for F-seq peak calling using `-l LENGTH`.
2. You could create wig files by setting `--wig`.
3. You could run with multiple threads using `-p THREAD`.

Example: 

```
call_peak.py rampage_TETSS_peak
```

### Step 3:

Calculate entropy for RAMPAGE peaks

```
Usage: entropy.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
```

* Input: the output folder created by `rm_pcr.py`
* Output: `rampage_entropy.txt` under the input folder

Format of `rampage_entropy.txt`:

The first thirteen columns of `rampage_entropy.txt` are the same as `rampage_peaks.txt`.

The additional three columns are listed below:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Entropy     | Entropy of RAMPAGE peak       |
| Positions_of_3'_ends | 3' ends of read pairs in the peak| 
|Read_counts_of_each_3'_end|Read counts of each 3' end  in the peak| 

Note:

1. You could run with multiple threads using `-p THREAD`.

Example: 

```
entropy.py rampage_TETSS_peak
```

### Step 4:

Assemble transcripts using RNA-seq (optional)

```
Usage: assemble.py [options] -g GTF <bam>...

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -g GTF --gtf=GTF               Gene annotation GTF file.
    -p THREAD --thread=THREAD      Threads. [default: 10]
    --dir=DIR                      Output directory. [default: ./]
    --prefix=PREFIX                Prefix for merged GTF.
                                   [default: merged_stringtie]
```
* Input:
    * BAM files of RNA-seq (`<bam>...`)
    * Gene annotation GTF file (`-g GTF`)
* Output: 
    * Assembled transcript GTF file for each RNA-seq dataset (`*stringtie.gtf`)
    * Merged transcript GTF file for all RNA-seq datasets (`merged_stringtie.gtf`)
    * Merged transcript BED file for all RNA-seq datasets (`merged_stringtie.bed`) 

Note:

1. This step is optional if there are no matched RNA-seq datasets available. If this step is skipped, the known gene annotation file could be used only in the downstream analysis.
2. You could run with multiple threads using `-p THREAD`.

Example: 
```
assemble.py -g gene_annotation.gtf rna_rep1.bam rna_rep2.bam 
```

### Step 5:

Fetch RAMPAGE peaks with TE

```
Usage: overlap_te.py [options] -t te <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -t te                          TE annotations (BED format).
    --extend length                TE extended length. [default: 50]
    --span span                    RAMPAGE span cutoff. [default: 1000]
    --entropy entropy              RAMPAGE entropy cutoff. [default: 2.5]
```

* Input: 
    * TE annotation BED file (`-t te`)
    * The output folder created by `entropy.py`
* Output:  `rampage_TE.txt` under the input folder
 

Format of `rampage_TE.txt`:

The first sixteen columns are the same as `rampage_entropy.txt`.

The additional seven columns are listed below:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome of TE              |
| Start       | Start of TE                   |
| End         | End of TE                     |
| Name        | Name of TE                    |
| Score       | 0                             |
| Strand      | Strand of TE                  |
| Direction   | Direction of TE in the peak   |

Note:

1. If using Repeatmask annotation file, you could download them from UCSC.
2. You could set the TE annotation extension length using `--extend length`.
3. You could set entropy cutoff using `--entropy entropy`.
4. You could set RAMPAGE effective length cutoff using `--span span`.


Example: 

```
overlap_te.py -t te.bed rampage_TETSS_peak
```

### Step 6:

Assign TE peaks to transcripts

```
Usage: assign_peak.py [options] -a ANNO -g ANNO <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -a ANNO                        Known gene annotation BED12 file.
    -g ANNO                        Assembled gene annotation BED12 file.
    -f fraction                    Cutoff of assignment fraction. [default: 0.5]
    -p THREAD                      Threads. [default: 10]
```

* Input: 
    * Known gene annotation BED file (`-a ANNO`)
    * Assembled gene annotation BED file (`-g ANNO`) 
    * The output folder created by `overlap_te.py`
* Output:  `rampage_TE_gene.txt` under the input folder

Format of `rampage_TE_gene.txt`:

The first twenty-three columns are the same as `rampage_TE.txt`.

The additional fifteen columns are listed below:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome of transcript      |
| Start       | Start of transcript           |
| End         | End of transcript             |
| Name        | Name of transcript            |
| Score       | 0                             |
| Strand      | Strand of transcript          |
| CDS start   | CDS start of transcript       |
| CDS end     | CDS end of transcript         |
| Score       | 0                             |
| Exon_counts | Exon counts of transcript     |
| Exon_sizes  | Exon sizes of transcript      |
| Exon_offsets| Exon offsets of transcript    |
| Annotation  | known or novel                |
| Distance    | Distance to gene annotaion    |
| Fraction    | Assignment fraction           |

Note:

1. If skipping the Step 4, you could use the known gene annotaion file for `-g ANNO` instead.
2. You could set assignment fraction cutoff using `-f fraction`.
3. You could run with multiple threads using `-p THREAD`.

Example: 

```
assign_peak.py -a gene_annotation.bed -g merged_stringtie.bed rampage_TETSS_peak
```

## License

Copyright (C) Xiao-Ou Zhang. See the [LICENSE](https://github.com/kepbod/rampage_te_tss/blob/master/LICENSE) file for license rights and limitations (MIT).