scAPA
================

This is a package and a shell script for alternative polyadenylation analysis of single-cell RNA-seq data. The shell script "scAPAscript.R" takes as input BAM files generated by tenx 3' RNA seq pipline and the results of cell clustering.

Prerequisites
=============

1.  Fasta and chromosome length files for human (hg19) or (and) mouse (mm10). Can be downloaded from [UCSC website](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/)
2.  [Samtools](http://www.htslib.org/download/) version 0.1.19-96b5f2294a or above.
3.  [Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) version v2.27.1-17-gf59480f or above.
4.  [Homer](http://homer.ucsd.edu/homer/introduction/install.html) version v4.8.3 or above. Uses homer's "MakeTagdirectory" and "findPeaks".
5.  [UMI\_tools](https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md)
6.  [Drop-seq\_tools-2.2.0](https://github.com/broadinstitute/Drop-seq/releases/tag/v2.2.0)
7.  R version 3.5.3, or above.
8.  Optional: [ChangePoint](https://sourceforge.net/projects/utr/files/), version 0.1.1

Installing
==========

-   Install all the above.
-   Install the R package devtools. In R:

<!-- -->

    ## install.packages("devtools")

-   Download scAPA.shell.script folder. Or use:

<!-- -->

    ## git clone https://github.com/ElkonLab/scAPA.git

-   Open configfile.txt. Fill in the paths to the fasta files, chromosome length files, and all the tools as instructed in the file. Leave "PATH" if the tool is in your PATH environment variable. Save your changes.

Using scAPA shell script
========================

Input Files
-----------

For an example of input files, see the example below. Put in one directory the following files (and only those files):

**1. Single-cell RNA seq BAM files**

BAM files generated by tenx 3' RNA seq pipline (aligned using cell ranger counts). Cell barcode is CB, and the molecular barcode is UB. The name of each BAM file should be of the following format: sample.name.bam. Avoid underscores in the names.

**2. Cluster annotations**

Tab-delimited two-column files with no header. The first column is the cell barcodes (CB). The second is the cluster assigned to the cell. Each file corresponds to a sample (BAM). Notice that the cell barcodes are the same as in the BAM. For example, if the end of the barcode in the BAM file is "-1", it should be the same in the text file. The names of the files should be of the following format: clusters\_sample.name.txt, where the sample names match those of the BAM files.

Usage
-----

Run the script as follow:

    ## path.to.R/bin/Rscript path.to.scAPAshellscript/scAPA.shell.script.R -p <path.to.files> -org <organism> -sp <path.to.script.dir> [options]

&lt;path.to.files&gt; is the path to the directory with BAM and cluster anotations tex files. Do not add / in the end of the path. -org The organism is either Mm for a mouse (mm10) or Hs for human (hg19). &lt;path.to.script.dir&gt; is the path to the scAPA.shell.script directory. For a list of options view the [Options.md](Options.md) file, or type:

    ## path.to.R/bin/Rscript path.to.scAPAshellscript/scAPA.shell.script.R --help

Example for usage
-----------------

The files used for this example are from [Lukassen et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132189/). The original fastq files can be obtined from [GSE104556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104556). **Please notice that the BAM files here have been down sampled, to contain 25% of the reads.** The result will differ from the results obtained from analyzing the full BAM files. Downlad the folder down.sampled.spermatogenesis from [this link](https://drive.google.com/open?id=1xK7lR2ECfJ-Cjb1f4bYnaA5JjYtdqrGA). After downloading and configuring scAPA.sell.script.R, use the following command:

    ## Rscript path.to.scAPAshellscript/scAPA.shell.script.R -p path/to/down.sampled.spermatogenesis -org Mm -sp path/to/scAPA.shell.script

By default, the script counts reads at a single-cell level. For faster analysis, count reads from cell clusters by using the "-sc false" option. The default number of cores to use is 30.

### Estimated Run Time

Running the script on the example time with 30 cores, took ttt. Bellow are estimated time (30 cores) for each stage of the analysis:

-   **Stage 1a** PCR duplicates removal - approximately 2 hours and 30 minutes.

-   **Stage 1b** Peak detection - aproximetly 45 minetes

-   **Stage 1c** Separating Peaks - approximately 3 hours.

-   **Stage 2** Quantifying the usage of each peak - approximately 4 hours.

-   **Stage 3** Peak filtering - approximately

-   **Stage 4** Statistical analysis

-   **Stage 5** Inferring global trends

Output Files
------------

For a full list of output files and their description, see the file [outputs.md](outputs.md) The following is a partial list of outputs:

-   **summary.UTR.txt (summary.UTR.txt)** A short summary of the 3'UTR peak analysis.

-   **ThreeUTR.peaks.txt (Intron.peaks.tt)** The 3'UTR peaks that passed filtering. The file contains the peak ID, gene symbol, ensemble ID, and genomic location.

-   **APA.events.txt** For 3'UTRs with more than one peak, gives the p-value, FDR corrected q-value of APA event tested across the clusters. Also given is the proximal PUI index for each cell cluster for each 3'UTR.

-   **UTRs.with.multiple.peaks.txt** a file containing the results of testing the differential usage of individual peaks across clusters. Only peaks from 3'UTR that came up significant in the analysis and had more than 2 peaks are analyzed. The file contains tables for each 3'UTR. The p-value and FDR q-value for each peak (from chi-square test for goodness of fit) are given. The peak usage index (PUI) for each cluster is given. Higher PUI in a cluster means higher usage of the peak in the cluster.

-   **Mean.Cell.PPUI.txt** gives the mean proximal PUI index for the cells analyzed.

Log files
---------

The script's logfile scAPA.script.log is created in scAPA directory. A directory Log.files is created containing the log files from the tools used (such as UMI tools).

Authors
=======

-   Eldad Shulman
-   Dr. Ran Elkon

License
=======

This project is licensed under the BSD 3 License - see the [LICENSE.md](LICENSE.md) file for details
