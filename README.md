**Plastid Genome Annotation**<br />
Copyright (C) 2017 Xiao-Jian Qu<br />

**Contact**<br />
quxiaojian@mail.kib.ac.cn<br />

**Prerequisites**<br />
BLAST+<br />
Perl<br />
Windows, Linux or Mac<br />

**General Introduction to PGA**<br />
PGA(Plastid Genome Annotation) is capable of annotating multiple plastid genomes using GenBank-format plastomes as reference. Three steps will be conducted to annotate plastomes: (1) extracting annotation features from GenBank-format reference plastomes, (2) blasting of annotation features against FASTA-format target plastomes, (3) generating GenBank-format files for each FASTA-format plastome sequence, and giving corresponding warning information for further manual check.<br />

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
        [-r -reference]    required: (default: reference) input directory name containing GenBank-format file(s) that from the same or close families.
        [-t -target]       required: (default: target) input directory name containing FASTA-format file(s) that you want to annotate.
        [-i -ir]           optional: (default: 1000) allowed minimum value for inverted-repeat (IR) length.
        [-d -degree]       optional: (default: 1) 1st (2nd, 3rd and so on) longest IR that you want to annotate.
        [-p -pidentity]    optional: (default: 40) any PCG with a TBLASTN percent identity less than this value will be listed in the log file and
                           will not be annotated.
        [-q -qcoverage]    optional: (default: 0.5,2) any PCG with a query coverage per annotated PCG less or greater than each of these two values (<1,>1)
                           will be listed in the log file.
        [-o -out]          optional: (default: gb) output directory name.
        [-f -form]         optional: (default: circular) circular or linear form for FASTA-format file.
        [-l -log]          optional: (default: warning) log file name containing warning information for annotated GenBank-format file(s).
```

**Test**<br />
(1) annotating your FASTA-format target plastomes.<br />
```
PGA.pl -r reference -t target
```
or
```
PGA.pl -r reference -t target -i 1000 -d 1 -p 40 -q 0.5,2 -o gb -f circular -l warning
```

(2) checking warning information in log file.<br />
(3) correcting your annotations according to warning information using Geneious.<br />

**Citation**<br />
If you use PGA in you scientific research, please cite:<br />
BLAST+<br />
Camacho C, et al. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10:421.<br />
PGA<br />
https://github.com/quxiaojian/PGA<br />
