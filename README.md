**Plastid Genome Annotation**<br />
Copyright (C) 2017 Xiao-Jian Qu<br />

**Contact**<br />
quxiaojian@mail.kib.ac.cn<br />

**Prerequisites**<br />
BLAST+<br />
Perl<br />
Windows, Linux or Mac<br />

**General Introduction to PGA**<br />
PGA(Plastid Genome Annotation) is capable of annotating multiple plastid genomes using genebank format plastomes as reference. Three steps will be conducted to annotate plastomes: (1) extracting annotation features from gb format reference plastomes, (2) blasting annotation features against fasta format plastome sequences, (3) generating gb format files for each fasta format sequence, and giving corresponding warning information need to be manually checked.<br />

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
        PGA.pl -r -s [-v -o -t -l]
        Copyright (C) 2017 Xiao-Jian Qu
        Please contact <quxiaojian@mail.kib.ac.cn>, if you have any bugs or questions.

        [-h -help]       help information.
        [-r -ref]        required: input directory name containing GenBank format file(s)
                         that from the same or close families. (default: reference)
        [-s -seq]        required: input directory name containing fasta format file(s)
                         that you want to annotate. (default: sequence)
        [-v -val]        optional: similarity value for BLAST results of amino acid. (default: 40)
        [-o -out]        optional: output directory name. (default: gb)
        [-t -type]       optional: circular or linear for fasta format file(s). (default: circular)
        [-l -log]        optional: log file name containing warning information
                         for annotated GenBank format file(s). (default: warning)
```

**Test**<br />
(1) annotating your fasta format plastome sequences.<br />
```
PGA.pl -r reference -s sequence
```
(2) checking warning information in log file.<br />
(3) correcting your annotations according to warning information using Geneious.<br />

**Citation**<br />
If you use PGA in you scientific research, please cite:<br />
BLAST+<br />
Camacho C, et al. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10:421.<br />
PGA<br />
https://github.com/quxiaojian/PGA<br />
