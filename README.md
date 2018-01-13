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

![PGA flowchart](https://github.com/quxiaojian/PGA/blob/master/PGA.png)

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

**Run Test**<br />
```
PGA.pl -r test/1/reference -t test/1/target
```
or
```
PGA.pl -r test/1/reference -t test/1/target -i 1000 -d 1 -p 40 -q 0.5,2 -o gb -f circular -l warning
```

**Input and Output**<br />
Annotation of the plastome of Rosa roxburghii using PGA. (a) ¡°Amborella_trichopoda.gb¡± shows the partial GenBank-formatted reference plastome of Amborella trichopoda, as revised from AJ506156. (b) ¡°Rosa_roxburghii.fasta¡± shows the partial FASTA-formatted target plastome of Rosa roxburghii, revised from NC_032038. (c) ¡°Rosa_roxburghii.gb¡± shows the output GenBank-formatted file containing partial annotation information for the target plastome of Rosa roxburghii. (d) ¡°warning.log¡± shows warning and statistical items during the annotation of the target plastome of Rosa roxburghii. The log file indicates the loss of the atpF intron in Rosa roxburghii. There are 114 total genes in the reference and target plastomes.<br />

```
(a) Amborella_trichopoda.gb
LOCUS       Amborella_trichopoda      162686 bp    DNA     circular UNA 08-JUN-2015
DEFINITION  Amborella trichopoda chloroplast genomic DNA, complete sequence.
ACCESSION   AJ506156
VERSION     AJ506156.2  GI:34481608
KEYWORDS    complete genome.
SOURCE      chloroplast Amborella trichopoda
  ORGANISM  Amborella trichopoda
            Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
            Spermatophyta; Magnoliophyta; basal Magnoliophyta; Amborellales;
            Amborellaceae; Amborella.
FEATURES          Location/Qualifiers
     source          1..162686
                    /organism="Amborella trichopoda"
                    /mol_type="genomic DNA"
     repeat_region    90951..117611
                    /note="inverted repeat region B; IRB repeat region"
                    /rpt_type="inverted"
     rRNA          complement(139284..142097)
                    /gene="rrn23"
                    /product="23S ribosomal RNA"
     gene           complement(139284..142097)
                    /gene="rrn23"
     tRNA          join(complement(4472..4508), complement(1840..1874))
                    /gene="trnK-UUU"
                    /product="tRNA-Lys"
     gene           complement(1840..4508)
                    /gene="trnK-UUU"
     CDS           join(complement(16186..16330), complement(14506..14915))
                    /gene="atpF"
                    /codon_start=1
                    /transl_table=11
                    /product="ATPase I subunit"
/translation="MKNVTDSFVSLGHWPSAGSFGFNTDIFATNPINLSVVLGVLIFF
                    GKGVLSDLLDNRKQRILSTIRNSEELRGGAIEQLEKARARLRKVEIEADEFRVNGYSE
                    IEREKSNLINAAYENLERLENYKNESIHFEQQRAMNQVRQRVFQQALQGALETLNSYL
                         NSELHLRTISANIGMLGTMKNITD"
     gene           complement(14506..16330)
                    /gene="atpF"
(b) Rosa_roxburghii.fasta
>Rosa_roxburghii
ATGGGCGAACGACGGGAATTGAACCCGCGCGTGGTGGATTCACAATCCACTGCCTTGATC
(c) Rosa_roxburghii.gb
LOCUS       Rosa_roxburghii  156749 bp    DNA     circular PLN 25-DEC-2017
FEATURES             Location/Qualifiers
     source          1..156749
                     /organism="Rosa_roxburghii"
                     /mol_type="genomic DNA"
     gene            106222..109027
                     /gene="rrn23"
     rRNA            106222..109027
                     /gene="rrn23"
                     /product="23S ribosomal RNA"
     gene            complement(1704..4278)
                     /gene="trnK-UUU"
     tRNA            join(complement(4242..4278), complement(1704..1738))
                     /gene="trnK-UUU"
                     /product="tRNA-Lys"
     gene            complement(12213..12767)
                     /gene="atpF"
     CDS             complement(12213..12767)
                     /gene="atpF"
                     /codon_start=1
                     /transl_table=11
                     /product="ATP synthase CF0 subunit I"
(d) warning.log
Rosa_roxburghii
Warning: atpF (negative one-intron PCG) lost intron!
Total number of genes in the reference plastome(s): 114.
Total number of genes annotated in the target plastome: 114.
All gene names from the reference plastome(s) that were not annotated in the target plastome:
```

**Recommendations for using PGA**<br />
(1) Users should carefully check the GenBank-formatted reference plastome. PGA is packaged with several properly annotated plastomes, and it is thus possible for users to use PGA to re-annotate a plastome that is intended to be used as a reference, in order to correct possible inaccuracies.<br />
(2) It is important that users select a reference plastome that contains sufficient numbers of annotated genes for the target taxa. The number of genes in the reference plastome(s) should equal or exceed the number in the target plastome(s). If the number of genes in the target is uncertain, it may be best to use multiple reference plastomes. The Amborella trichopoda (AJ506156) and Zamia furfuracea (JX416857) plastomes included within PGA are examples of plastomes that contain the highest gene numbers among known angiosperms and gymnosperms, and as such it is recommended that they be included as references during PGA runs.<br />
(3) We do not recommend annotating highly incomplete plastomes using a complete reference plastome, because BLAST may annotate some genes redundantly (i.e., BLAST may return hits for genes that were not sequenced or are otherwise absent in the incomplete plastome, resulting in spurious annotations). To annotate highly incomplete plastomes or plastome segments, we recommend using progressiveMauve (as implemented in Mauve 2.4.0; Darling et al., 2010) to align the incomplete plastome to the reference plastome, followed by the use of the corresponding homologous block of the reference plastome as the reference for annotation in PGA.<br />
(4) We suggest that users carefully check highly divergent or otherwise unusual target plastomes for incorrect annotations. This is particularly important for plastomes with a high degree of gene loss, pseudogenization or sequence divergence.<br />

**Boundary Detection Algorithms**<br />
Three algorithms are applied to (1) determine start and stop codons, (2) locate intron-exon boundaries and detect intron loss, and (3) identify the boundaries of the Inverted Repeat (IR). Following two figures show the first two algorithms, respectively. IR boundary annotation is accomplished via a self-BLASTN search. Two parameters can be adjusted to determine the IR boundaries: minimum allowed IR length (default = 1000) and the first (second, third, and so on) longest inverted repeat that will be annotated as the IR (default = 1).<br />
![GBDA](https://github.com/quxiaojian/PGA/blob/master/GBDA.png)

![IBDA](https://github.com/quxiaojian/PGA/blob/master/IBDA.png)

**Citation**<br />
If you use PGA in you scientific research, please cite:<br />
BLAST+<br />
Camacho C, et al. 2009. BLAST+: architecture and applications. BMC Bioinformatics 10:421.<br />
PGA<br />
https://github.com/quxiaojian/PGA<br />
