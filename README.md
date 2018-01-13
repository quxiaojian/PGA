**Plastid Genome Annotation**<br />
Copyright (C) 2018 Xiao-Jian Qu<br />

**Contact**<br />
quxiaojian@mail.kib.ac.cn<br />

**Prerequisites**<br />
BLAST+<br />
Perl<br />
Windows, Linux or Mac<br />

**General Introduction to PGA**<br />
PGA (Plastid Genome Annotation), a standalone command line tool, can perform rapid, accurate, and flexible batch annotation of newly generated target plastomes based on well-annotated reference plastomes. In contrast to current existing tools, PGA uses reference plastomes as the query and unannotated target plastomes as the subject to locate genes, which we refer to as the reverse query-subject BLAST search approach. PGA accurately identifies gene and intron boundaries as well as intron loss. The program outputs GenBank-formatted files as well as a log file to assist users in verifying annotations.<br />

Following six steps will be conducted to annotate plastomes: (1) Preparation of GenBank-formatted reference plastomes; (2) Preparation of FASTA-formatted target plastomes; (3) Reference database generation; (4) BLAST search; (5) Determining feature boundaries; (6) Generating GenBank and log files.<br />

![PGA flowchart](https://github.com/quxiaojian/PGA/blob/master/PGA.tiff)

**Preparations**<br />

(1) download BLAST+ software [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download), and put it in PATH.<br />
```
vim ~/.bashrc
export PATH=/home/xxx/blast-2.7.0+/bin:$PATH
source ~/.bashrc
```
(2) download this repository to your local computer, and put it in PATH. Make it read, write and executable.<br />
```
git clone https://github.com/quxiaojian/PGA.git
vim ~/.bashrc
export PATH=/home/xxx/PGA:$PATH
source ~/.bashrc
chmod a+rwx PGA.pl
```

You can test PGA.pl by type PGA.pl, which will show the usage information.<br />
```
Usage:
    PGA.pl -r -t [-i -d -p -q -o -f -l]
    Copyright (C) 2018 Xiao-Jian Qu
    Please contact <quxiaojian@mail.kib.ac.cn>, if you have any bugs or questions.

    [-h -help]         help information.
    [-r -reference]    required: (default: reference) input directory name containing GenBank-formatted file(s) that from the same or close families.
    [-t -target]       required: (default: target) input directory name containing FASTA-formatted file(s) that will be annotated.
    [-i -ir]           optional: (default: 1000) minimum allowed inverted-repeat (IR) length.
    [-d -degree]       optional: (default: 1) the first (second, third, and so on) longest inverted repeat that will be annotated as the IR.
    [-p -pidentity]    optional: (default: 40) any PCGs with a TBLASTN percent identity less than this value will be listed in the log file and
                       will not be annotated.
    [-q -qcoverage]    optional: (default: 0.5,2) any PCGs with a query coverage per annotated PCG less or greater than each of these two values (<1,>1)
                       will be listed in the log file.
    [-o -out]          optional: (default: gb) output directory name.
    [-f -form]         optional: (default: circular) circular or linear form for FASTA-formatted file.
    [-l -log]          optional: (default: warning) log file name containing warning information for annotated GenBank-formatted file(s).
```

**Test**<br />
(1) annotating your FASTA-formatted target plastomes.<br />
```
PGA.pl -r reference -t target
```
or
```
PGA.pl -r reference -t target -i 1000 -d 1 -p 40 -q 0.5,2 -o gb -f circular -l warning
```

(2) checking warning information in log file.<br />
(3) correcting your annotations using Geneious according to warning information in log file.<br />

**Boundary Detection Algorithms**<br />
Three algorithms are applied to (1) determine start and stop codons, (2) locate intron-exon boundaries and detect intron loss, and (3) identify the boundaries of the Inverted Repeat (IR). Following two figures show the first two algorithms, respectively. IR boundary annotation is accomplished via a self-BLASTN search. Two parameters can be adjusted to determine the IR boundaries: minimum allowed IR length (default = 1000) and the first (second, third, and so on) longest inverted repeat that will be annotated as the IR (default = 1).<br />
![GBDA](https://github.com/quxiaojian/PGA/blob/master/GBDA.tiff)

![IBDA](https://github.com/quxiaojian/PGA/blob/master/IBDA.tiff)

**Citation**<br />
If you use PGA in you scientific research, please cite:<br />
BLAST+<br />
Camacho C, et al. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10:421.<br />
PGA<br />
https://github.com/quxiaojian/PGA<br />
