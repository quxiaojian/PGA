##Plastid Genome Annotation<br />
Copyright (C) 2017 Xiao-Jian Qu<br />

##Contact<br />
quxiaojian@mail.kib.ac.cn<br />
Notes: This script is beta version, some aspects need to be improved.<br />

##Prerequisites<br />
GNU<br />
Perl<br />
Linux and Windows<br />

##General Introduction to PGA<br />
PGA(Plastid Genome Annotation) is capable of annotating multiple plastid genomes using published genebank format files as reference. Three steps will be conducted to annotate plastome: (1) extracting annotation information from gb format references, (2) blasting annotation information to fasta format sequences, (3) generating gb format files for fasta format sequences and giving warning information need to be manually checked.<br />

##Preparations<br />

You can test PGA.pl by type ~/PATH/TO/PGA.pl, which will show the usage information.<br />
```
Usage:
    PGA.pl -r -s [-o -t -l]
    Copyright (C) 2017 Xiao-Jian Qu
    Please contact <quxiaojian@mail.kib.ac.cn>, if you have any bugs or questions.

    [-h -help]       help information.
    [-r -ref]        required: input directory name containing gb reference file(s) that from the same or close family,order etc.(default: reference)
    [-s -seq]        required: input directory name containing fasta sequence files(s) that you want to annotate.(default: sequence)
    [-o -out]        optional: output directory name containing annotated genebank(gb) file(s).(default: gb)
    [-t -type]       optional: circular or linear for fasta sequence files(s).(default: circular)
    [-l -log]        optional: log file name containing warning information for annotated genebank(gb) file(s).(default: warning)
```

##Tutorial<br />
**First**, annotating your fasta format plastome sequences.<br />
```
PGA.pl -r reference -s sequence
```
**Second**, checking your annotations for each sequence using Geneious.<br />
