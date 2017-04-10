#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Find;
use Time::HiRes qw(time);
use Data::Dumper;
no warnings "uninitialized";
$|=1;

my $global_options=&argument();
my $reference_directory=&default("reference","ref");
my $sequence_directory=&default("sequence","seq");
my $type=&default("circular","type");
my $output_directory=&default("gb","out");
my $log=&default("warning","log");

my $now1=time;
my $now2=&gettime;
print "$now2 || Begin extracting annotations from reference!\n";

print "\nPGA.pl Plastid Genome Annotation
Copyright (C) 2016 Xiao-Jian Qu
Email: quxiaojian\@mail.kib.ac.cn\n\n";

system("del/f/s/q $output_directory") if (-e $output_directory);
mkdir ($output_directory) if (!-e $output_directory);

############################################################
## extract_gene_coding_CDS
############################################################
my $pattern=".gb";
my @filenames;
find(\&target1,$reference_directory);
sub target1{
    if (/$pattern/){
        push @filenames,"$File::Find::name";
    }
    return;
}

while (@filenames) {
	my $filename_gb=shift @filenames;
   	my $filename_base=$filename_gb;
	$filename_base=~ s/(.*).gb/$1/g;
	my $latin_name=substr ($filename_base,rindex($filename_base,"\/")+1);
	open(my $in_gb,"<",$filename_gb);
	open(my $out_gb,">","$filename_base\_temp1");
	while (<$in_gb>){
		$_=~ s/\r\n/\n/g;
		if ($_=~ /\),\n/){
			$_=~ s/\),\n/\),/g;
		}elsif($_=~ /,\n/){
			$_=~ s/,\n/,/g;
		}
		print $out_gb $_;
	}
	close $in_gb;
	close $out_gb;

	open(my $in_gbk,"<","$filename_base\_temp1");
	open(my $out_gbk,">","$filename_base\_temp2");
	while (<$in_gbk>){
		$_=~ s/,\s+/,/g;
		print $out_gbk $_;
	}
	close $in_gbk;
	close $out_gbk;


	############################################################
	## extract_bed_gene
	############################################################
	if (!-e "$filename_base\_gene.bed"){
		#generate_bed_file
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		open (my $in_genbank,"<","$filename_base\_temp2");
		while (<$in_genbank>){
			chomp;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$latin_name;
				$length=$row_array[2];
			}elsif(/ {5}gene {12}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				$element=$species_name.":".$1.":".$element;
				push @genearray,$element;
				$element=();
				$mark=0;
			}
		}
		close $in_genbank;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();

		#put_fasta_sequence_in_array
	    my $flag=0;
	    my @sequence;
		my (@fas1,@fas2);
		open(my $in_genebank,"<","$filename_base\_temp2");
	    while (<$in_genebank>){
	        if ($_=~ /ORIGIN/){
	            $flag=1;
	        }
	        if ($_=~ /\/\//){
	            $flag=2;
	        }
	        if ($flag==1){
	            next if ($_=~ /ORIGIN/);
	            push @sequence,$_;
	        }
	    }
		close $in_genebank;
		foreach (@sequence){
			chomp;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
	    my $fas1=join "",@fas1;
	    my (@fasta1,@fasta2);
	    push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			chomp;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9];
		}
		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif (defined $STRAND2{$_} ne "") {# rps12_gene,including rps12+1_gene and rps12+2_gene
				if (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif (($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}
			}
		}
		@output1=();

		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();

		#output_bed_file
		open (my $out_bed,">","$filename_base\_gene.bed");
		foreach (@output3){
			print $out_bed $_;
		}
		close $out_bed;

		#extract_gene
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$filename_base\_gene.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		foreach (@output3){
			chomp;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			if ($STRAND4{$cnt2} eq "-") {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				print $out_coding ">".$GENE2{$cnt2}."_gene_".$SP2{$cnt2}."\n".$rev_com."\n";
			}elsif($STRAND4{$cnt2} eq "+"){
				print $out_coding ">".$GENE2{$cnt2}."_gene_".$SP2{$cnt2}."\n".$str."\n";
			}
		}
		close $out_coding;

		#generate_IGS_ranges
		my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last0,$last1,$last2,@output4);
		my $cnt3=0;
		foreach (@output3){
			chomp;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START5{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
	        $last0=$_-1;
			$last1=$_+1;
			$last2=$_+2;
			next if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) < ($END5{$last1}-1)));
			next if (($_ > 1) and (($END5{$_}+1) < ($END5{$last0}-1)) and (($END5{$_}+1) < ($START5{$last2}-1)));
			if ((($END5{$_}+1) >= ($START5{$last1}-1)) and (($END5{$_}+1) >= ($END5{$last1}-1))){
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND5{$_}."/".$STRAND5{$last2}."\t".($END5{$_}+1)."\t".($START5{$last2}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last2}."\n";
	        }else{
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND5{$_}."/".$STRAND5{$last1}."\t".($END5{$_}+1)."\t".($START5{$last1}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last1}."\n";
	        }
		}
		foreach (keys %SP3){
			if ($_==$cnt3){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();

		#extract_IGS
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$filename_base\_IGS.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		foreach (@output4){
			chomp;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			print $out_noncoding ">".$GENE4{$cnt4}."_gene_".$SP4{$cnt4}."\n".$str."\n";
		}
		@output4=();
		close $out_noncoding;
	    unlink "$filename_base\_IGS.fasta";
    }


	############################################################
	## extract_bed_coding
	############################################################
	if (!-e "$filename_base\_coding.bed"){
		#generate_bed_file
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		open (my $in_genbank,"<","$filename_base\_temp2");
		while (<$in_genbank>){
			chomp;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$latin_name;
				$length=$row_array[2];
			}elsif(/ {5}CDS {13}/ or / {5}tRNA {12}/ or / {5}rRNA {12}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				$element=$species_name.":".$1.":".$element;
				push @genearray,$element;
				$element=();
				$mark=0;
			}
		}
		close $in_genbank;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();

		#put_fasta_sequence_in_array
	    my $flag=0;
	    my @sequence;
		my (@fas1,@fas2);
		open(my $in_genebank,"<","$filename_base\_temp2");
	    while (<$in_genebank>){
	        if ($_=~ /ORIGIN/){
	            $flag=1;
	        }
	        if ($_=~ /\/\//){
	            $flag=2;
	        }
	        if ($flag==1){
	            next if ($_=~ /ORIGIN/);
	            push @sequence,$_;
	        }
	    }
		close $in_genebank;
		foreach (@sequence){
			chomp;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
	    my $fas1=join "",@fas1;
	    my (@fasta1,@fasta2);
	    push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			chomp;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1},$STRAND3{$cnt1},$START3{$cnt1},$END3{$cnt1},$TYPE3{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
		}
		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "")) {
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "")) {# rps12_coding, inlcuding rps12+1_coding, rps12+2-1_coding and rps12+2-2_coding
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}
			}
		}
		@output1=();

		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();

		#output_bed_file
		open (my $out_bed,">","$filename_base\_coding.bed");
		foreach (@output3){
			print $out_bed $_;
		}
		close $out_bed;

		#extract_coding_regions
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$filename_base\_coding.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}
		foreach (@output3){
			chomp;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
			if ($STRAND4{$cnt2} eq "-") {
				my $rev_com=reverse $str;
				$rev_com=~ tr/ACGTacgt/TGCAtgca/;
				print $out_coding ">".$GENE2{$cnt2}."_coding_".$SP2{$cnt2}."\n".$rev_com."\n";
			}elsif($STRAND4{$cnt2} eq "+"){
				print $out_coding ">".$GENE2{$cnt2}."_coding_".$SP2{$cnt2}."\n".$str."\n";
			}
		}
		close $out_coding;

		#generate_noncoding_ranges
		my (%SP3,%GENE3,%STRAND5,%START5,%END5,%TYPE5,$last,@output4);
		my $cnt3=0;
		foreach (@output3){
			chomp;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND5{$cnt3},$START5{$cnt3},$END5{$cnt3},$TYPE5{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START5{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND5{$_}."\t"."1"."\t".($START5{$_}-1)."\t"."?"."/".$TYPE5{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
			$last=$_+1;
			next if ($END5{$_}+1) >= ($START5{$last}-1);
			push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last}."\t".$STRAND5{$_}."/".$STRAND5{$last}."\t".($END5{$_}+1)."\t".($START5{$last}-1)."\t".$TYPE5{$_}."/".$TYPE5{$last}."\n";
		}
		foreach (keys %SP3){
			if ($_==$cnt3){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND5{$_}."/"."?"."\t".($END5{$_}+1)."\t".$length."\t".$TYPE5{$_}."/"."?"."\n";
			}
		}
		@output3=();

		#extract_noncoding_regions
		my (%SP4,%GENE4,%STRAND6,%START6,%END6,%TYPE6,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$filename_base\_noncoding.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		foreach (@output4){
			chomp;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND6{$cnt4},$START6{$cnt4},$END6{$cnt4},$TYPE6{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START6{$cnt4}-1),($END6{$cnt4}-$START6{$cnt4}+1));
			print $out_noncoding ">".$GENE4{$cnt4}."_coding_".$SP4{$cnt4}."\n".$str."\n";
		}
		@output4=();
		close $out_noncoding;
	    unlink "$filename_base\_noncoding.fasta";
    }


	############################################################
	## extract_bed_CDS
	############################################################
	if (!-e "$filename_base\_CDS.bed"){
		#generate_bed_file
		my (@row_array,$species_name,$length,$element,@genearray,@output1);
		my $mark=0;
		open (my $in_genbank,"<","$filename_base\_temp2");
		while (<$in_genbank>){
			chomp;
			@row_array=split /\s+/,$_;
			if (/^LOCUS/i){
				$species_name=$latin_name;
				$length=$row_array[2];
			}elsif(/ {5}CDS {13}/){
				if ($row_array[2]=~ /^\d+..\d+$/){
					$row_array[2]="\+\t$row_array[2]\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~/^complement\((\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^complement\(join\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\((\d+..\d+),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),complement\((\d+..\d+)\),complement\((\d+..\d+)\)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t-\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^order\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\(complement\((\d+..\d+)\),(\d+..\d+),(\d+..\d+)\)$/){
					$row_array[2]="-\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t+\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}elsif($row_array[2]=~ /^join\((\d+..\d+),(\d+..\d+),complement\((\d+..\d+)\)\)$/) {
					$row_array[2]="+\t$1\t$row_array[1]\t+\t$2\t$row_array[1]\t-\t$3\t$row_array[1]";
					$row_array[2]=~ s/\../\t/g;
				}
				$element=$row_array[2];
				$mark=1;
			}elsif(/^\s+\/gene="(.*)"/ and $mark == 1){
				$element=$species_name.":".$1.":".$element;
				push @genearray,$element;
				$element=();
				$mark=0;
			}
		}
		close $in_genbank;
		foreach (@genearray){
			my @array=split /:/,$_;
			push @output1,"$array[0]\t$array[1]\t$array[2]\n";
		}
		@row_array=();
		@genearray=();

		#put_fasta_sequence_in_array
	    my $flag=0;
	    my @sequence;
		my (@fas1,@fas2);
		open(my $in_genebank,"<","$filename_base\_temp2");
	    while (<$in_genebank>){
	        if ($_=~ /ORIGIN/){
	            $flag=1;
	        }
	        if ($_=~ /\/\//){
	            $flag=2;
	        }
	        if ($flag==1){
	            next if ($_=~ /ORIGIN/);
	            push @sequence,$_;
	        }
	    }
		close $in_genebank;
		foreach (@sequence){
			chomp;
			$_=~ s/\s*//g;
			$_=~ s/\d+//g;
			push @fas1,$_;
		}
	    my $fas1=join "",@fas1;
	    my (@fasta1,@fasta2);
	    push @fasta1,$species_name,$fas1;
		@fasta2=@fasta1;

		#edit_bed_file
		my (%SP1,%GENE1,%STRAND1,%START1,%END1,%TYPE1,%STRAND2,%START2,%END2,%TYPE2,%STRAND3,%START3,%END3,%TYPE3,@output2);
		my $cnt1=0;
		foreach (@output1) {
			chomp;
			$cnt1++;
			($SP1{$cnt1},$GENE1{$cnt1},$STRAND1{$cnt1},$START1{$cnt1},$END1{$cnt1},$TYPE1{$cnt1},$STRAND2{$cnt1},$START2{$cnt1},$END2{$cnt1},$TYPE2{$cnt1},$STRAND3{$cnt1},$START3{$cnt1},$END3{$cnt1},$TYPE3{$cnt1})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
		}
		foreach (1..$cnt1) {
			if (defined $STRAND2{$_} eq "") {
				push @output2,($SP1{$_}."\t".$GENE1{$_}."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} eq "")) {
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_})){
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
				}
			}elsif ((defined $STRAND2{$_} ne "") and (defined $STRAND3{$_} ne "")) {# rps12_CDS, including rps12+1_CDS, rps12+2-1_CDS and rps12+2-2_CDS
				if (($STRAND1{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "-") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($START1{$_} > $START2{$_}) and ($START2{$_} > $START3{$_}) and ($GENE1{$_} ne "rps12")){# non-rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-3"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START2{$_}) and ($START2{$_} < $START3{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "-") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}elsif(($STRAND1{$_} eq "+") and ($STRAND2{$_} eq "+") and ($START1{$_} < $START3{$_}) and ($START3{$_} < $START2{$_}) and ($GENE1{$_} eq "rps12")){# rps12
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+1"."\t".$STRAND1{$_}."\t".$START1{$_}."\t".$END1{$_}."\t".$TYPE1{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-2"."\t".$STRAND2{$_}."\t".$START2{$_}."\t".$END2{$_}."\t".$TYPE2{$_});
					push @output2,($SP1{$_}."\t".$GENE1{$_}."+2-1"."\t".$STRAND3{$_}."\t".$START3{$_}."\t".$END3{$_}."\t".$TYPE3{$_});
				}
			}
		}

		#sort_bed_file
		my $col=3;
		my (%sort,@output3);
		foreach (@output2){
			my @row=split /\t/,$_;
			$sort{$_}=$row[$col];
		}
		foreach (sort {$sort{$a} <=> $sort{$b}} keys %sort){
			push @output3,"$_\n";
		}
		@output2=();

		#output_bed_file
		open (my $out_bed,">","$filename_base\_CDS.bed");
		foreach (@output3){
			print $out_bed $_;
		}
		close $out_bed;

		#extract_CDS
		my (%SP2,%GENE2,%STRAND4,%START4,%END4,%TYPE4,%STRAND5,%START5,%END5,%TYPE5,%STRAND6,%START6,%END6,%TYPE6,$seq1);
		my $cnt2=0;
		open(my $out_coding,">","$filename_base\_CDS.fasta");
		while (@fasta1){
			my $header=shift @fasta1;
			$seq1=shift @fasta1;
		}

		foreach (@output1){
			chomp;
			$cnt2++;
			($SP2{$cnt2},$GENE2{$cnt2},$STRAND4{$cnt2},$START4{$cnt2},$END4{$cnt2},$TYPE4{$cnt2},$STRAND5{$cnt2},$START5{$cnt2},$END5{$cnt2},$TYPE5{$cnt2},$STRAND6{$cnt2},$START6{$cnt2},$END6{$cnt2},$TYPE6{$cnt2})=(split /\s+/,$_)[0,1,2,3,4,5,6,7,8,9,10,11,12,13];
			if (defined $STRAND5{$cnt2} eq "") {
	        	my $str1=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	            if ($STRAND4{$cnt2} eq "-") {
	                my $rev_com1=reverse $str1;
	                $rev_com1=~ tr/ACGTacgt/TGCAtgca/;
	                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com1."\n";
	            }elsif($STRAND4{$cnt2} eq "+"){
	                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str1."\n";
	            }
	#        }elsif((defined $STRAND5{$cnt2} ne "") and (defined $STRAND6{$cnt2} eq "")) {
	#            if (($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} < $START5{$cnt2})) {
	#                my $str2=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
	#                my $rev_com2=reverse $str2;
	#                $rev_com2=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com2."\n";
	#            }elsif(($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} > $START5{$cnt2})) {
	#                my $str3=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $rev_com3=reverse $str3;
	#                $rev_com3=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com3."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} < $START5{$cnt2})){
	#                my $str4=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str4."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} > $START5{$cnt2})){
	#                my $str5=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str5."\n";
	#            }
	#        }elsif ((defined $STRAND5{$cnt2} ne "") and (defined $STRAND6{$cnt2} ne "")) {# generating rps12+1_CDS and rps12+2_CDS
	#            if (($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2}) and ($GENE2{$cnt2} ne "rps12")) {# non-rps12
	#                my $str6=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
	#                my $rev_com4=reverse $str6;
	#                $rev_com4=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com4."\n";
	#            }elsif(($STRAND4{$cnt2} eq "-") and ($START4{$cnt2} > $START5{$cnt2}) and ($START5{$cnt2} > $START6{$cnt2}) and ($GENE2{$cnt2} ne "rps12")) {# non-rps12
	#                my $str7=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $rev_com5=reverse $str7;
	#                $rev_com5=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com5."\n";
	#            }elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "+") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str8=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str9=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
	#                my $rev_com6=reverse $str8;
	#                $rev_com6=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."+1_CDS_".$SP2{$cnt2}."\n".$rev_com6."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."+2_CDS_".$SP2{$cnt2}."\n".$str9."\n";
	#            }elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "-") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str10=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str11=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
	#                my $rev_com7=reverse $str10;
	#                $rev_com7=~ tr/ACGTacgt/TGCAtgca/;
	#                my $rev_com8=reverse $str11;
	#                $rev_com8=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."+1_CDS_".$SP2{$cnt2}."\n".$rev_com7."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."+2_CDS_".$SP2{$cnt2}."\n".$rev_com8."\n";
	#            }elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "+") and ($START4{$cnt2} < $START6{$cnt2}) and ($START6{$cnt2} < $START5{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str12=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str13=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
	#                my $rev_com9=reverse $str12;
	#                $rev_com9=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."+1_CDS_".$SP2{$cnt2}."\n".$rev_com9."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."+2_CDS_".$SP2{$cnt2}."\n".$str13."\n";
	#            }elsif(($STRAND4{$cnt2} eq "-") and ($STRAND5{$cnt2} eq "-") and ($START4{$cnt2} < $START6{$cnt2}) and ($START6{$cnt2} < $START5{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str14=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str15=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
	#                my $rev_com10=reverse $str14;
	#                $rev_com10=~ tr/ACGTacgt/TGCAtgca/;
	#                my $rev_com11=reverse $str15;
	#                $rev_com11=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."+1_CDS_".$SP2{$cnt2}."\n".$rev_com10."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."+2_CDS_".$SP2{$cnt2}."\n".$rev_com11."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2}) and ($GENE2{$cnt2} ne "rps12")){# non-rps12
	#                my $str16=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str16."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($START4{$cnt2} > $START5{$cnt2}) and ($START5{$cnt2} > $START6{$cnt2}) and ($GENE2{$cnt2} ne "rps12")){# non-rps12
	#                my $str17=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str17."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "-") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str18=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str19=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
	#                my $rev_com12=reverse $str19;
	#                $rev_com12=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."+1_CDS_".$SP2{$cnt2}."\n".$str18."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."+2_CDS_".$SP2{$cnt2}."\n".$rev_com12."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "+") and ($START4{$cnt2} < $START5{$cnt2}) and ($START5{$cnt2} < $START6{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str20=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str21=substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1)).substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1));
	#                print $out_coding ">".$GENE2{$cnt2}."+1_CDS_".$SP2{$cnt2}."\n".$str20."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."+2_CDS_".$SP2{$cnt2}."\n".$str21."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "-") and ($START4{$cnt2} < $START6{$cnt2}) and ($START6{$cnt2} < $START5{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str22=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str23=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
	#                my $rev_com13=reverse $str23;
	#                $rev_com13=~ tr/ACGTacgt/TGCAtgca/;
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str22."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$rev_com13."\n";
	#            }elsif(($STRAND4{$cnt2} eq "+") and ($STRAND5{$cnt2} eq "+") and ($START4{$cnt2} < $START6{$cnt2}) and ($START6{$cnt2} < $START5{$cnt2}) and ($GENE2{$cnt2} eq "rps12")) {# rps12
	#				my $str24=substr($seq1,($START4{$cnt2}-1),($END4{$cnt2}-$START4{$cnt2}+1));
	#                my $str25=substr($seq1,($START6{$cnt2}-1),($END6{$cnt2}-$START6{$cnt2}+1)).substr($seq1,($START5{$cnt2}-1),($END5{$cnt2}-$START5{$cnt2}+1));
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str24."\n";
	#                print $out_coding ">".$GENE2{$cnt2}."_CDS_".$SP2{$cnt2}."\n".$str25."\n";
	#            }
	        }
		}
		close $out_coding;
	    @output1=();

		#generate_intergenic_ranges
		my (%SP3,%GENE3,%STRAND7,%START7,%END7,%TYPE7,$last0,$last1,$last2,@output4);
		my $cnt3=0;
		foreach (@output3){
			chomp;
			$cnt3++;
			($SP3{$cnt3},$GENE3{$cnt3},$STRAND7{$cnt3},$START7{$cnt3},$END7{$cnt3},$TYPE7{$cnt3})=(split /\s+/,$_)[0,1,2,3,4,5];
		}
		foreach (keys %SP3){
			if ($_==1 and $START7{$_}!=1){
				unshift @output4,$SP3{$_}."\t"."start".'-'.$GENE3{$_}."\t"."?"."/".$STRAND7{$_}."\t"."1"."\t".($START7{$_}-1)."\t"."?"."/".$TYPE7{$_}."\n";
			}
		}
		foreach (1..($cnt3-1)) {
	        $last0=$_-1;
			$last1=$_+1;
			$last2=$_+2;
			next if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) < ($END7{$last1}-1)));
			next if (($_ > 1) and (($END7{$_}+1) < ($END7{$last0}-1)) and (($END7{$_}+1) < ($START7{$last2}-1)));
			if ((($END7{$_}+1) >= ($START7{$last1}-1)) and (($END7{$_}+1) >= ($END7{$last1}-1))){
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last2}."\t".$STRAND7{$_}."/".$STRAND7{$last2}."\t".($END7{$_}+1)."\t".($START7{$last2}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last2}."\n";
	        }else{
	    		push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'.$GENE3{$last1}."\t".$STRAND7{$_}."/".$STRAND7{$last1}."\t".($END7{$_}+1)."\t".($START7{$last1}-1)."\t".$TYPE7{$_}."/".$TYPE7{$last1}."\n";
	        }
		}
		foreach (keys %SP3){
			if ($_==$cnt3){
				push @output4,$SP3{$_}."\t".$GENE3{$_}.'-'."end"."\t".$STRAND7{$_}."/"."?"."\t".($END7{$_}+1)."\t".$length."\t".$TYPE7{$_}."/"."?"."\n";
			}
		}
		@output3=();

		#extract_intergenic
		my (%SP4,%GENE4,%STRAND8,%START8,%END8,%TYPE8,$seq2);
		my $cnt4=0;
		open(my $out_noncoding,">","$filename_base\_intergenic.fasta");
		while (@fasta2){
			my $header=shift @fasta2;
			$seq2=shift @fasta2;
		}
		foreach (@output4){
			chomp;
			$cnt4++;
			($SP4{$cnt4},$GENE4{$cnt4},$STRAND8{$cnt4},$START8{$cnt4},$END8{$cnt4},$TYPE8{$cnt4})=(split /\s+/,$_)[0,1,2,3,4,5];
			my $str=substr($seq2,($START8{$cnt4}-1),($END8{$cnt4}-$START8{$cnt4}+1));
			print $out_noncoding ">".$GENE4{$cnt4}."_CDS_".$SP4{$cnt4}."\n".$str."\n";
		}
		@output4=();
		close $out_noncoding;
	    unlink "$filename_base\_intergenic.fasta";
    }
	unlink "$filename_base\_temp1";
	unlink "$filename_base\_temp2";
}


############################################################
## combine_gene_coding_CDS
############################################################
my $dir;
opendir($dir,$reference_directory);
my @input_directory=readdir $dir;
close $dir;
unlink ("gene.fasta");
unlink ("coding.fasta");
unlink ("CDS.fasta");

foreach my $file (@input_directory){
	if($file=~/_gene.fasta/){
		open (my $input_gene,"<","$reference_directory/$file");
		open (my $output_gene,">>","gene.fasta");
		while(<$input_gene>){
			print $output_gene $_;
		}
		close $input_gene;
		close $output_gene;
	}
	if($file=~/_coding.fasta/){
		open (my $input_coding,"<","$reference_directory/$file");
		open (my $output_coding,">>","coding.fasta");
		while(<$input_coding>){
			print $output_coding $_;
		}
		close $input_coding;
		close $output_coding;
	}
	if($file=~/_CDS.fasta/){
		open (my $input_CDS,"<","$reference_directory/$file");
		open (my $output_CDS,">>","CDS.fasta");
		while(<$input_CDS>){
			print $output_CDS $_;
		}
		close $input_CDS;
		close $output_CDS;
	}
}

open (my $in_gene,"<","gene.fasta");
open (my $in_coding,"<","coding.fasta");
open (my $in_CDS,"<","CDS.fasta");
open (my $out_all1,">","all1.fasta");
while (<$in_gene>){
	print $out_all1 $_;
}
while (<$in_coding>){
	print $out_all1 $_;
}
while (<$in_CDS>){
	print $out_all1 $_;
}
close $in_gene;
close $in_coding;
close $in_CDS;
close $out_all1;

#from A-UGC to trnA-UGC
open (my $in_all1,"<","all1.fasta");
open (my $out_all2,">","all2.fasta");
my ($header_in_all1,$sequence_in_all1);
while (defined ($header_in_all1=<$in_all1>) && defined ($sequence_in_all1=<$in_all1>)){
	if ($header_in_all1=~ />(.|..)-/){
		$header_in_all1=~ s/>/>trn/g;
		print $out_all2 $header_in_all1.$sequence_in_all1;
	}else{
		print $out_all2 $header_in_all1.$sequence_in_all1;
	}
}
close $in_all1;
close $out_all2;
unlink("all1.fasta");

#from CDS to CDS_aa(no intron PCG)
#from coding to coding_aa(intron PCG)
open (my $in_all2,"<","all2.fasta");
open (my $out_all3,">","all3.fasta");
my %hash_in_all2=("---"=>"-","TAA"=>"*","TAG"=>"*","TGA"=>"*","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S","TTC"=>"F","TTT"=>"F","TTA"=>"L","TTG"=>"L","TAC"=>"Y","TAT"=>"Y","TGC"=>"C","TGT"=>"C","TGG"=>"W","CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P","CAC"=>"H","CAT"=>"H","CAA"=>"Q","CAG"=>"Q","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R","ATA"=>"I","ATC"=>"I","ATT"=>"I","ATG"=>"M","ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T","AAC"=>"N","AAT"=>"N","AAA"=>"K","AAG"=>"K","AGC"=>"S","AGT"=>"S","AGA"=>"R","AGG"=>"R","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E","GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G");
my ($header_in_all2,$sequence_in_all2);
while (defined ($header_in_all2=<$in_all2>) && defined ($sequence_in_all2=<$in_all2>)){
	chomp ($header_in_all2,$sequence_in_all2);
	my $length=length $sequence_in_all2;
	my $aa;

	if ($header_in_all2=~ /_CDS/){
		for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
		#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
			my $codon=substr ($sequence_in_all2,$i,3);
			$codon=uc $codon;
			if (exists $hash_in_all2{$codon}){
				$aa.=$hash_in_all2{$codon};
			}else{
				$aa.="X";
				my $j=$i+1;
				#print STDERR "Bad codon $codon in position $j of species $header_in_all2!\n";
			}
		}
		$header_in_all2=~ s/_CDS/_CDS_aa/g;
		print $out_all3 "$header_in_all2\n$aa\n";
	}elsif ($header_in_all2=~ /^>(?!trn)(.+)-(\d)_coding/){
		if ($2 == 1) {
			for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
			#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
				my $codon=substr ($sequence_in_all2,$i,3);
				$codon=uc $codon;
				if (exists $hash_in_all2{$codon}){
					$aa.=$hash_in_all2{$codon};
				}else{
					$aa.="X";
					my $j=$i+1;
					#print STDERR "Bad codon $codon in position $j of species $header_in_all2!\n";
				}
			}
		}elsif ($2 == 2) {
			if ($length % 3 == 0) {
				for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
				#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
					my $codon=substr ($sequence_in_all2,$i,3);
					$codon=uc $codon;
					if (exists $hash_in_all2{$codon}){
						$aa.=$hash_in_all2{$codon};
					}else{
						$aa.="X";
						my $j=$i+1;
						#print STDERR "Bad codon $codon in position $j of species $header_in_all2!\n";
					}
				}
			}elsif ($length % 3 == 1) {
				for (my $i=1;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
				#for (my $i=1;$i<($length-3);$i+=3){# delete stop codon
					my $codon=substr ($sequence_in_all2,$i,3);
					$codon=uc $codon;
					if (exists $hash_in_all2{$codon}){
						$aa.=$hash_in_all2{$codon};
					}else{
						$aa.="X";
						my $j=$i+1;
						#print STDERR "Bad codon $codon in position $j of species $header_in_all2!\n";
					}
				}
			}elsif ($length % 3 == 2) {
				for (my $i=2;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
				#for (my $i=2;$i<($length-3);$i+=3){# delete stop codon
					my $codon=substr ($sequence_in_all2,$i,3);
					$codon=uc $codon;
					if (exists $hash_in_all2{$codon}){
						$aa.=$hash_in_all2{$codon};
					}else{
						$aa.="X";
						my $j=$i+1;
						#print STDERR "Bad codon $codon in position $j of species $header_in_all2!\n";
					}
				}
			}
		}elsif ($2 == 3) {
			for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
			#for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
				my $codon=substr ($sequence_in_all2,$i,3);
				$codon=uc $codon;
				if (exists $hash_in_all2{$codon}){
					$aa.=$hash_in_all2{$codon};
				}else{
					$aa.="X";
					my $j=$i+1;
					#print STDERR "Bad codon $codon in position $j of species $header_in_all2!\n";
				}
			}
		}
		$header_in_all2=~ s/_coding/_coding_aa/g;
		print $out_all3 "$header_in_all2\n$aa\n";
	}else{
		print $out_all3 "$header_in_all2\n$sequence_in_all2\n";
	}
}
close $in_all2;
close $out_all3;

#delete xxx_coding of no intron gene(including PCG and RNA):
#if exists accD_coding and accD_gene,delete accD_coding(no intron PCG);
#if exists trnQ-UUG_coding and trnQ-UUG_gene,delete trnQ-UUG_coding(no intron tRNA);
#if exists rrn16_coding and rrn16_gene,delete rrn16_coding(no intron rRNA);
open (my $in_all3,"<","all3.fasta");
open (my $out_reference1,">","reference1.fasta");
open (my $out_reference2,">","reference2.fasta");
open (my $out_reference3,">","reference3.fasta");
open (my $out_reference4,">","reference4.fasta");

my ($header_in_all3,$sequence_in_all3,%hash1_in_all3,%hash2_in_all3,%hash3_in_all3,%hash4_in_all3);
while (defined ($header_in_all3=<$in_all3>) and defined ($sequence_in_all3=<$in_all3>)){
	chomp ($header_in_all3,$sequence_in_all3);
	if ($header_in_all3=~ /^>(trn)(.+)-1_coding/ or $header_in_all3=~ /^>(trn)(.+)-2_coding/ or $header_in_all3=~ /^>(trn)(.+)-3_coding/){
		$hash1_in_all3{$header_in_all3}=$sequence_in_all3;
	}
	if ($header_in_all3=~ /^>(trn)(.+)_gene/ or $header_in_all3=~ /^>(rrn)(.+)_gene/){
		$hash1_in_all3{$header_in_all3}=$sequence_in_all3;
	}
	if ($header_in_all3=~ /^>(?!trn)(.+)_gene/ and $header_in_all3=~ /^>(?!rrn)(.+)_gene/){
		$hash2_in_all3{$header_in_all3}=$sequence_in_all3;
	}
	if ($header_in_all3=~ /^>(.+)_CDS_aa/){
		$hash3_in_all3{$header_in_all3}=$sequence_in_all3;
	}
	if ($header_in_all3=~ /^>(.+)_coding_aa/){
		$hash4_in_all3{$header_in_all3}=$sequence_in_all3;
	}
}
foreach my $h1 (sort keys %hash1_in_all3){
	print $out_reference1 "$h1\n$hash1_in_all3{$h1}\n";
}
foreach my $h2 (sort keys %hash2_in_all3){
	print $out_reference2 "$h2\n$hash2_in_all3{$h2}\n";
}
foreach my $h3 (sort keys %hash3_in_all3){
	print $out_reference3 "$h3\n$hash3_in_all3{$h3}\n";
}
foreach my $h4 (sort keys %hash4_in_all3){
	print $out_reference4 "$h4\n$hash4_in_all3{$h4}\n";
}
close $in_all3;
close $out_reference1;
close $out_reference2;
close $out_reference3;
close $out_reference4;
unlink("all3.fasta");
unlink ("gene.fasta");
unlink ("coding.fasta");
unlink ("CDS.fasta");
%hash1_in_all3=();
%hash2_in_all3=();
%hash3_in_all3=();
%hash4_in_all3=();

#short_exon_of_protein-coding-gene_or_tRNA
open (my $input_ref,"<","all2.fasta");
my ($h_ref,$s_ref,%hash_exon,%hash_exon_length,%hash_exon_sequence);
while (defined ($h_ref=<$input_ref>) and defined ($s_ref=<$input_ref>)){
	chomp ($h_ref,$s_ref);
	my $len=length $s_ref;
	if ($h_ref=~ /^>((?!trn)(.+)-1_coding)/){
		$hash_exon{$1}{$len}++;
		$hash_exon_length{$1}=$len;
		$hash_exon_sequence{$1}=$s_ref;
	}elsif($h_ref=~ /^>((?!trn)(.+)-2_coding)/){
		$hash_exon{$1}{$len}++;
		$hash_exon_length{$1}=$len;
		$hash_exon_sequence{$1}=$s_ref;
	}elsif($h_ref=~ /^>((?!trn)(.+)-3_coding)/){
		$hash_exon{$1}{$len}++;
		$hash_exon_length{$1}=$len;
		$hash_exon_sequence{$1}=$s_ref;
	}
}
close $input_ref;
unlink("all2.fasta");
#print Dumper \%hash_exon;
#print Dumper \%hash_exon_length;
#print Dumper \%hash_exon_sequence;
my $now3=&gettime;
print "$now3 || Finish extracting annotations from reference!\n";




############################################################
## sequence_dealing_in_batch_mode
############################################################
my $pattern1=".fasta";
my $pattern2=".fas";
my $pattern3=".fa";
my @sequence_filenames;
find(\&target2,$sequence_directory);
sub target2{
    if (/$pattern1/ or /$pattern2/ or /$pattern3/){
        push @sequence_filenames,"$File::Find::name";
    }
    return;
}

my $j=0;
while (@sequence_filenames) {
	$j++;
	my $input_fasta=shift @sequence_filenames;
	my $output_fasta=substr($input_fasta,0,index ($input_fasta,"\."));
	open(my $input_ag,"<",$input_fasta);
	open(my $output_ag,">",$output_fasta);
	my $row_ag=<$input_ag>;
	print $output_ag $row_ag;
	while ($row_ag=<$input_ag>){
		chomp $row_ag;
		if ($row_ag=~ /^>/) {
			print $output_ag "\n".$row_ag."\n";
		}else{
			print $output_ag $row_ag;
		}
	}
	print $output_ag "\n";
	close $input_ag;
	close $output_ag;

	#fasta_sequence
	open (my $in_fasta,"<",$output_fasta);
	my @fasta;
	while (<$in_fasta>){
		chomp;
		push @fasta,$_;
	}
	close $in_fasta;
	my ($header,$sequence,$length);
	while (@fasta){
		my $head=shift @fasta;
		$head=~ s/ +/_/g;
		$head=~ s/_$//g;
		$header=$1 if ($head=~ /^>(.+)$/);
		$sequence=shift @fasta;
		$length=length $sequence;
	}


	############################################################
	## blast_reference_to_sequence
	############################################################
	#rRNA and no intron tRNA(_gene)
	#intron tRNA(_gene and -1/-2_coding)
	#no intron PCG (_gene and _CDS_aa)
	#intron PCG(_gene and -1/-2/-3_coding_aa)(short exon of rpl16,petB and petD)

	system ("makeblastdb.exe -in $output_fasta -hash_index -dbtype nucl");

	#-max_hsps 1(nucleotide,RNA,including _gene RNA and -1/-2_coding tRNA)
	system ("blastn.exe -task blastn -query reference1.fasta -db $output_fasta -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out blast_reference1");# -evalue 0.001 or 0.01
	#-max_hsps 1(nucleotide,PCG,including _gene PCG)
	system ("blastn.exe -task blastn -query reference2.fasta -db $output_fasta -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out blast_reference2");# -evalue 0.001 or 0.01
	#-max_hsps 1(amino acid,no intron PCG,including _CDS_aa PCG)
	system ("tblastn.exe -task tblastn -query reference3.fasta -db $output_fasta -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out blast_reference3");# -evalue 0.001 or 0.01
	#-max_hsps 1(amino acid,intron PCG,including -1/-2/-3_coding_aa PCG)
	system ("tblastn.exe -task tblastn -query reference4.fasta -db $output_fasta -outfmt 6 -max_hsps 1 -max_target_seqs 1 -out blast_reference4");# -evalue 0.001 or 0.01

	unlink("$output_fasta");
	unlink ("reference_temp");

	open (my $in_bt_ref1,"<","blast_reference1");
	open (my $in_bt_ref2,"<","blast_reference2");
	open (my $in_bt_ref3,"<","blast_reference3");
	open (my $in_bt_ref4,"<","blast_reference4");
	open (my $out_bt_ref,">>","reference_temp");
	while(<$in_bt_ref1>){
		print $out_bt_ref $_;
	}
	while(<$in_bt_ref2>){
		print $out_bt_ref $_;
	}
	while(<$in_bt_ref3>){
		print $out_bt_ref $_;
	}
	while(<$in_bt_ref4>){
		print $out_bt_ref $_;
	}
	close $in_bt_ref1;
	close $in_bt_ref2;
	close $in_bt_ref3;
	close $in_bt_ref4;
	close $out_bt_ref;

	open (my $in_ref,"<","reference_temp");
	open (my $out_ref,">","reference.tab");
	my @array_reference;
	while (<$in_ref>){
		chomp;
		my ($qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore)=split (/\s+/,$_);
		push @array_reference,[$qseqid,$sseqid,$pident,$length,$mismatch,$gapopen,$qstart,$qend,$sstart,$send,$evalue,$bitscore];
	}

	foreach my $item (sort {$a->[0] cmp $b->[0] or $a->[1] cmp $b->[1] or $b->[2] <=> $a->[2] or $b->[3] <=> $a->[3]} @array_reference){
		print $out_ref "$item->[0]\t$item->[1]\t$item->[2]\t$item->[3]\t$item->[4]\t$item->[5]\t$item->[6]\t$item->[7]\t$item->[8]\t$item->[9]\t$item->[10]\t$item->[11]\n";
	}
	close $in_ref;
	close $out_ref;
	unlink ("blast_reference1");
	unlink ("blast_reference2");
	unlink ("blast_reference3");
	unlink ("blast_reference4");
	unlink ("reference_temp");


	############################################################
	## generate_annotation_table
	############################################################
	open (my $input_reference,"<","reference.tab");
	open (my $output_annotation_tab,">","$header.tab");
	my ($row_reference,%hash_tab,%hash1_tab,%hash2_tab,%hash3_tab,%hash4_tab,%hash_gene);
	my $i=0;
	while (defined ($row_reference=<$input_reference>)){
		my @array=split /\t/,$row_reference;
		my $contig=$array[1];
		%hash1_tab=();
		if ($array[0]=~ /((.+)_CDS_aa)/){
			$hash_gene{$1}=1;
			$hash1_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash1_tab)==1);
		}
		%hash2_tab=();
		if ($array[0]=~ /((.+)_gene)/){
			$hash_gene{$1}=1;
			$hash2_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash2_tab)==1);
		}
		%hash3_tab=();
		if ($array[0]=~ /((.+)-(\d)_coding_aa)/){
			$hash_gene{$1}=1;
			$hash3_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash3_tab)==1);
		}
		%hash4_tab=();
		if ($array[0]=~ /((.+)-(\d)_coding)(?!_aa)/){
			$hash_gene{$1}=1;
			$hash4_tab{$2}++;
			$hash_tab{$contig}{$2}{$1}{$array[8]."\t".$array[9]}++ if ((keys %hash4_tab)==1);
		}
	}
	close $input_reference;
	#print Dumper \%hash_tab;

	my @gene_coding_CDS;
	foreach (keys %hash_gene){
		push @gene_coding_CDS,$_;
	}
	foreach my $contig_name (sort keys %hash_tab){
		print $output_annotation_tab "$contig_name\n";
		foreach my $name (sort keys %{$hash_tab{$contig_name}}){
			print $output_annotation_tab "$name\t";
			my $cnt1=0;
			foreach my $gene_coding_CDS (sort {$b cmp $a} keys %{$hash_tab{$contig_name}{$name}}){
				print $output_annotation_tab "$gene_coding_CDS\t" if ($cnt1==0);
				print $output_annotation_tab "\t\t$gene_coding_CDS\t" if ($cnt1>0);
				my @keys;
				foreach my $key1 (sort {$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$b} <=> $hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$a}} keys %{$hash_tab{$contig_name}{$name}{$gene_coding_CDS}}){
					push @keys,$key1;
				}
				my %large_hash;
				if ($gene_coding_CDS=~ /((.+)_CDS_aa)/){
					my $gene1=$2."-1_coding";
					my $gene2=$2."-2_coding";
					my $gene3=$2."-3_coding";
					if ((grep {$gene1 ne $_} @gene_coding_CDS) and (grep {$gene2 ne $_} @gene_coding_CDS)){
						$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
					}
					if(((grep {$gene1 eq $_} @gene_coding_CDS) or (grep {$gene2 eq $_} @gene_coding_CDS)) and (grep {$gene3 ne $_} @gene_coding_CDS)){
						$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
						if (defined $keys[1]){
							$large_hash{$keys[1]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[1]};
						}
					}
					if(grep {$gene3 eq $_} @gene_coding_CDS){
						$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
						if (defined $keys[1]){
							$large_hash{$keys[1]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[1]};
						}
						if (defined $keys[2]){
							$large_hash{$keys[2]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[2]};
						}
					}
				}else{
					$large_hash{$keys[0]}=$hash_tab{$contig_name}{$name}{$gene_coding_CDS}{$keys[0]};
				}
				#print Dumper \%large_hash;
				foreach my $key6 (sort keys %large_hash){
						print $output_annotation_tab "$key6\t$large_hash{$key6}\t";
				}
				$cnt1++;
				print $output_annotation_tab "\n";
				%large_hash=();
			}
		}
	}
	close $output_annotation_tab;
	unlink("reference.tab");
	%hash_tab=();
	%hash1_tab=();
	%hash2_tab=();
	%hash3_tab=();
	%hash4_tab=();
	%hash_gene=();




	############################################################
	## generate_gb_file_from_annotation_table
	############################################################
	#codon_table
	my %hash_codon=("---"=>"-","TAA"=>"*","TAG"=>"*","TGA"=>"*","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S","TTC"=>"F","TTT"=>"F","TTA"=>"L","TTG"=>"L","TAC"=>"Y","TAT"=>"Y","TGC"=>"C","TGT"=>"C","TGG"=>"W","CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","CCA"=>"P","CCC"=>"P","CCG"=>"P","CCT"=>"P","CAC"=>"H","CAT"=>"H","CAA"=>"Q","CAG"=>"Q","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R","ATA"=>"I","ATC"=>"I","ATT"=>"I","ATG"=>"M","ACA"=>"T","ACC"=>"T","ACG"=>"T","ACT"=>"T","AAC"=>"N","AAT"=>"N","AAA"=>"K","AAG"=>"K","AGC"=>"S","AGT"=>"S","AGA"=>"R","AGG"=>"R","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GAC"=>"D","GAT"=>"D","GAA"=>"E","GAG"=>"E","GGA"=>"G","GGC"=>"G","GGG"=>"G","GGT"=>"G");

	#product_of_protein-coding-gene_and_tRNA,e.g., $hash_tRNA{trnQ-UUG}=tRNA-Gln
	open (my $in_product,"<","product.txt");
	my %hash_product;
	while (<$in_product>){
		chomp;
		my @product=split /\t/,$_;
		$hash_product{$product[0]}=$product[1];
	}
	close $in_product;
	#print Dumper \%hash_product;

	open (my $in_annotation,"<","$header.tab");
	my $contig=<$in_annotation>;
	chomp $contig;
	my %hash;
	while (<$in_annotation>){
		my @array=split /\t/,$_;
		if ($_!~ /^\t/){
			my $gene=$array[0];
			if ($array[1]=~ /((.+)_gene)/){
				$hash{$gene}{$1}{$array[2]."\t".$array[3]}=$array[4];
			}
		}elsif($_=~ /^\t/){
			if ($array[2]=~ /((.+)_CDS_aa)/){
				$hash{$2}{$1}{$array[3]."\t".$array[4]}=$array[5];
				#$hash{$2}{$1}{$array[6]."\t".$array[7]}=$array[8] if (defined $array[8]);
			}
			if ($array[2]=~ /((.+)-(\d)_coding_aa)/){
				$hash{$2}{$1}{$array[3]."\t".$array[4]}=$array[5];
			}
			if ($array[2]=~ /((.+)-(\d)_coding(?!_aa))/){
				$hash{$2}{$1}{$array[3]."\t".$array[4]}=$array[5];
			}
		}
	}
	close $in_annotation;
	unlink("$header.tab");
	#print Dumper \%hash;

	open (my $out_annotation,">","$output_directory/$header.gb");
	open (my $logfile,">>","$output_directory/$log.log");
	my $time=&getdate;
	print $logfile "$header\n";
	print $out_annotation "LOCUS       $header  $length bp    DNA     $type PLN $time"."\n";
	print $out_annotation "FEATURES             Location/Qualifiers"."\n";
	print $out_annotation "     source          "."1..$length"."\n";
	print $out_annotation "                     /organism=\"$header\""."\n";
	print $out_annotation "                     /mol_type=\"genomic DNA\""."\n";

	foreach my $name (sort keys %hash){
		############################################################
		## no-intron RNA
		############################################################
		if ((keys %{$hash{$name}} == 1)){
			foreach my $gene_coding_CDS (keys %{$hash{$name}}){
				my ($start,$end);
				foreach my $position (keys %{$hash{$name}{$gene_coding_CDS}}){
					($start,$end)=(split /\t/,$position)[0,1];
					print $out_annotation "     "."gene"."            ".$start."..".$end."\n" if ($start < $end);
					print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n" if ($start > $end);
					print $out_annotation "                     "."/gene=\"$name\""."\n";
				}
				if (($gene_coding_CDS=~ /(((.+)-(\D+))_gene)/) or ($gene_coding_CDS=~ /trnR_gene/) or ($gene_coding_CDS=~ /trnA_gene/)){# tRNA
					print $out_annotation "     "."tRNA"."            ".$start."..".$end."\n" if ($start < $end);
					print $out_annotation "     "."tRNA"."            "."complement(".$end."..".$start.")\n" if ($start > $end);
					print $out_annotation "                     "."/gene=\"$2\""."\n";
					print $out_annotation "                     "."/product=\"".$hash_product{$2}."\""."\n";
				}
				if($gene_coding_CDS=~ /((rrn(\d+\.?\d*))_gene)/){# rRNA
					print $out_annotation "     "."rRNA"."            ".$start."..".$end."\n" if ($start < $end);
					print $out_annotation "     "."rRNA"."            "."complement(".$end."..".$start.")\n" if ($start > $end);
					print $out_annotation "                     "."/gene=\"$2\""."\n";
					print $out_annotation "                     "."/product=\"".$hash_product{$2}."\""."\n";
				}
				if(($gene_coding_CDS=~ /((.+)_gene)/) and $2!~ /rrn/ and $2!~ /trn/ and $2!~ /rps12/){# pseudogene (e.g., ycf15,ycf68 etal.)
					print $out_annotation "     "."CDS"."             ".$start."..".$end."\n" if ($start < $end);
					print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n" if ($start > $end);
					print $out_annotation "                     "."/gene=\"$2\""."\n";
					print $out_annotation "                     "."/codon_start=1"."\n";
					print $out_annotation "                     "."/transl_table=11"."\n";
					print $out_annotation "                     "."/product=\"".$hash_product{$2}."\""."\n";
					#print $out_annotation "                     "."/translation=\"$2\""."\n";
					print $logfile "Warning: $name does not have CDS type in reference and is most likely pseudogene!\n";
				}
			}
		}

		############################################################
		## no-intron protein-coding-gene
		############################################################
		if((keys %{$hash{$name}} == 2)){
			my (@position1,@position2);
			foreach my $gene_coding_CDS (sort {$b cmp $a} keys %{$hash{$name}}){
				if ($gene_coding_CDS=~ /((.+)_gene)/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position1,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /((.+)_CDS_aa)/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position2,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}
			}
			while (@position1 and @position2){
				my $position1=shift @position1;
				my $position2=shift @position2;
				my ($start1,$end1,$number1)=(split /\t/,$position1)[0,1,2];# _gene
				my ($start2,$end2,$number2)=(split /\t/,$position2)[0,1,2];# _CDS_aa
				if (($start1 == $start2) and ($end1 == $end2)){# identical PCG boundary for _gene and _CDS_aa
					if ($start1 < $end1){# positive
						my $str=substr($sequence,($start1-1),($end1-$start1+1));
						my $length=length $str;
						my $forward_start_codon=substr($str,0,3);
						my $forward_stop_codon=substr($str,($length-3),3);
						my $aa;
						#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
							my $codon=substr ($str,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa.=$hash_codon{$codon};
							}else{
								$aa.="X";
								my $j=$i+1;
								print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa=split //,$aa;
						if (($length % 3==0) and (!grep {$_=~ /\*/} @aa) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA")) or ($name eq "rps12+1")){# standard start and stop codon
							print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."CDS"."             ".$start1."..".$end1."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/codon_start=1"."\n";
							print $out_annotation "                     "."/transl_table=11"."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $out_annotation "                     "."/translation=\"$aa\""."\n";
						}elsif((grep {$_=~ /\*/} @aa) or (($forward_start_codon ne "ATG") or ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") or ($forward_stop_codon ne "TAG") or ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon
							my $gene_length=($end1-$start1+1);
							my ($start_left,$start_right,$end_left,$end_right);
							if ($start2>=9000) {
								$start_left=$start2-9000;
							}elsif ($start2<9000){
								$start_left=$start2-(int($start2/3)-1)*3;
							}
							$start_right=$start2+60;
							$end_left=$end2-int($gene_length/3)*3;
							if (($length-$end2)>=9000) {
								$end_right=$end2+9000;
							}elsif (($length-$end2)<9000){
								$end_right=$end2+(int(($length-$end2)/3)-1)*3;
							}

							my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
							my $str2=substr($sequence,$end_left,($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $aa1;#left range for start codon
							#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($str1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
									push @star1,$_;
								}
								if ($aa1[$_]=~ /M/){
									push @start1,$_;
								}
							}
							my $left_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
							}
							my $left_start_position=$start2[0];
							my $aa2;#right range for stop codon
							#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($str2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $right_star_position=shift @star2;
							if ((defined $left_start_position) and (defined $right_star_position)){
								my $start=$start_left+$left_start_position * 3;
								my $end=$end_left+($right_star_position * 3+3);
								my $str=substr($sequence,($start-1),($end-$start+1));
								my $length=length $str;
								my $aa;
								#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($str,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
							}elsif ((!defined $left_start_position) and (defined $right_star_position)){
								my $end=$end_left+($right_star_position * 3+3);
								print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             ".$start2."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
								print $logfile "Warning: $name (positive no-intron PCG) has alternative start codon!\n";
							}
						}
					}elsif($start1 > $end1){# negative
						my $str=substr($sequence,($end1-1),($start1-$end1+1));
						my $rev_com=reverse $str;
						$rev_com=~ tr/ACGTacgt/TGCAtgca/;
						my $length=length $str;
						my $reverse_start_codon=substr($rev_com,0,3);
						my $reverse_stop_codon=substr($rev_com,($length-3),3);
						my $aa;
						#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
							my $codon=substr ($rev_com,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa.=$hash_codon{$codon};
							}else{
								$aa.="X";
								my $j=$i+1;
								print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa=split //,$aa;
						if (($length % 3==0) and (!grep {$_=~ /\*/} @aa) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA")) or ($name eq "rps12+1")){# standard start and stop codon
							print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1.")\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/codon_start=1"."\n";
							print $out_annotation "                     "."/transl_table=11"."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $out_annotation "                     "."/translation=\"$aa\""."\n";
						}elsif((grep {$_=~ /\*/} @aa) or (($reverse_start_codon ne "ATG") or ($reverse_start_codon ne "GTG")) and (($reverse_stop_codon ne "TAA") or ($reverse_stop_codon ne "TAG") or ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon
							my $gene_length=($start1-$end1+1);
							my ($start_left,$start_right,$end_left,$end_right);
							$start_left=$start2-60;
							if (($length-$start2)>=9000) {
								$start_right=$start2+9000;
							}elsif (($length-$start2)<9000){
								$start_right=$start2+(int(($length-$start2)/3)-1)*3;
							}
							if ($end2>=9000) {
								$end_left=$end2-9000;
							}elsif ($end2<9000){
								$end_left=$end2-(int($end2/3)-1)*3;
							}
							$end_right=$end2+int($gene_length/3)*3;

							my $str1=substr($sequence,($start_left),($start_right-$start_left));
							my $str2=substr($sequence,($end_left-1),($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $aa1;#left range for start codon
							#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
									push @star1,$_;
								}
								if ($aa1[$_]=~ /M/){
									push @start1,$_;
								}
							}
							my $right_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
							}
							my $right_start_position=$start2[0];
							my $aa2;#right range for stop codon
							#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $left_star_position=shift @star2;
							if ((defined $right_start_position) and (defined $left_star_position)){
								my $start=$start_right-$right_start_position * 3;
								my $end=$end_right-($left_star_position * 3+3);
								my $str=substr($sequence,($end-1),($start-$end+1));
								my $length=length $str;
								my $rev_com=reverse $str;
								$rev_com=~ tr/ACGTacgt/TGCAtgca/;
								my $aa;
								#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
							}elsif ((!defined $right_start_position) and (defined $left_star_position)){
								my $end=$end_right-($left_star_position * 3+3);
								print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start2.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
								print $logfile "Warning: $name (negative no-intron PCG) has alternative start codon!\n";
							}
						}
					}
				}elsif(($start1 != $start2) or ($end1 != $end2)){# non-identical PCG boundary for _gene and _CDS_aa
					if ($start2 < $end2){# positive
						my $str1=substr($sequence,($start1-1),($end1-$start1+1));
						my $length1=length $str1;
						my $forward_start_codon1=substr($str1,0,3);
						my $forward_stop_codon1=substr($str1,($length1-3),3);
						my $str2=substr($sequence,($start2-1),($end2-$start2+1));
						my $length2=length $str2;
						my $forward_start_codon2=substr($str2,0,3);
						my $forward_stop_codon2=substr($str2,($length2-3),3);
						my $aa1;
						#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
							my $codon=substr ($str1,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa1.=$hash_codon{$codon};
							}else{
								$aa1.="X";
								my $j=$i+1;
								print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa1=split //,$aa1;
						my $aa2;
						#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
							my $codon=substr ($str2,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa2.=$hash_codon{$codon};
							}else{
								$aa2.="X";
								my $j=$i+1;
								print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa2=split //,$aa2;
						if (($length2 % 3==0) and (!grep {$_=~ /\*/} @aa2) and (($forward_start_codon2 eq "ATG") or ($forward_start_codon2 eq "GTG")) and (($forward_stop_codon2 eq "TAA") or ($forward_stop_codon2 eq "TAG") or ($forward_stop_codon2 eq "TGA"))){# standard start and stop codon for _CDS_aa
							print $out_annotation "     "."gene"."            ".$start2."..".$end2."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."CDS"."             ".$start2."..".$end2."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/codon_start=1"."\n";
							print $out_annotation "                     "."/transl_table=11"."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
						}elsif(($length1 % 3==0) and (!grep {$_=~ /\*/} @aa1) and (($forward_start_codon1 eq "ATG") or ($forward_start_codon1 eq "GTG")) and (($forward_stop_codon1 eq "TAA") or ($forward_stop_codon1 eq "TAG") or ($forward_stop_codon1 eq "TGA"))){# standard start and stop codon for _gene
							print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."CDS"."             ".$start1."..".$end1."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/codon_start=1"."\n";
							print $out_annotation "                     "."/transl_table=11"."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
						}elsif((grep {$_=~ /\*/} @aa2) or (($forward_start_codon2 ne "ATG") and ($forward_start_codon2 ne "GTG")) or (($forward_stop_codon2 ne "TAA") and ($forward_stop_codon2 ne "TAG") and ($forward_stop_codon2 ne "TGA"))){# non-standard start or stop codon for _CDS_aa
							my $gene_length=($end2-$start2+1);
							my ($start_left,$start_right,$end_left,$end_right);
							if ($start2>=9000) {
								$start_left=$start2-9000;
							}elsif ($start2<9000){
								$start_left=$start2-(int($start2/3)-1)*3;
							}
							$start_right=$start2+60;
							$end_left=$end2-int($gene_length/3)*3;
							if (($length-$end2)>=9000) {
								$end_right=$end2+9000;
							}elsif (($length-$end2)<9000){
								$end_right=$end2+(int(($length-$end2)/3)-1)*3;
							}

							my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
							my $str2=substr($sequence,($end_left),($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $aa1;#left range for start codon
							#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($str1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
									push @star1,$_;
								}
								if ($aa1[$_]=~ /M/){
									push @start1,$_;
								}
							}
							my $left_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
							}
							my $left_start_position=$start2[0];
							my $aa2;#right range for stop codon
							#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($str2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $right_star_position=shift @star2;
							if ((defined $left_start_position) and (defined $right_star_position)){
								my $start=$start_left+$left_start_position * 3;
								my $end=$end_left+($right_star_position * 3+3);
								my $str=substr($sequence,($start-1),($end-$start+1));
								my $length=length $str;
								my $aa;
								#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($str,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             ".$start."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
							}elsif ((!defined $left_start_position) and (defined $right_star_position)){
								my $end=$end_left+($right_star_position * 3+3);
								print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             ".$start2."..".$end."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
								print $logfile "Warning: $name (positive no-intron PCG) has alternative start codon!\n";
							}
						}
					}elsif($start2 > $end2){# negative
						my $str1=substr($sequence,($end1-1),($start1-$end1+1));
						my $rev_com1=reverse $str1;
						$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
						my $length1=length $str1;
						my $reverse_start_codon1=substr($rev_com1,0,3);
						my $reverse_stop_codon1=substr($rev_com1,($length1-3),3);
						my $str2=substr($sequence,($end2-1),($start2-$end2+1));
						my $rev_com2=reverse $str2;
						$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
						my $length2=length $str2;
						my $reverse_start_codon2=substr($rev_com2,0,3);
						my $reverse_stop_codon2=substr($rev_com2,($length2-3),3);
						my $aa1;
						#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
							my $codon=substr ($rev_com1,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa1.=$hash_codon{$codon};
							}else{
								$aa1.="X";
								my $j=$i+1;
								print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa1=split //,$aa1;
						my $aa2;
						#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
						for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
							my $codon=substr ($rev_com2,$i,3);
							$codon=uc $codon;
							if (exists $hash_codon{$codon}){
								$aa2.=$hash_codon{$codon};
							}else{
								$aa2.="X";
								my $j=$i+1;
								print "Bad codon $codon in position $j of gene $name in species $contig!\n";
							}
						}
						my @aa2=split //,$aa2;
						if (($length2 % 3==0) and (!grep {$_=~ /\*/} @aa2) and (($reverse_start_codon2 eq "ATG") or ($reverse_start_codon2 eq "GTG")) and (($reverse_stop_codon2 eq "TAA") or ($reverse_stop_codon2 eq "TAG") or ($reverse_stop_codon2 eq "TGA"))){# standard start and stop codon for _CDS_aa
							print $out_annotation "     "."gene"."            "."complement(".$end2."..".$start2.")\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."CDS"."             "."complement(".$end2."..".$start2.")\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/codon_start=1"."\n";
							print $out_annotation "                     "."/transl_table=11"."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $out_annotation "                     "."/translation=\"$aa2\""."\n";
						}elsif(($length1 % 3==0) and (!grep {$_=~ /\*/} @aa1) and (($reverse_start_codon1 eq "ATG") or ($reverse_start_codon1 eq "GTG")) and (($reverse_stop_codon1 eq "TAA") or ($reverse_stop_codon1 eq "TAG") or ($reverse_stop_codon1 eq "TGA"))){# standard start and stop codon for _gene
							print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."CDS"."             "."complement(".$end1."..".$start1.")\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/codon_start=1"."\n";
							print $out_annotation "                     "."/transl_table=11"."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $out_annotation "                     "."/translation=\"$aa1\""."\n";
						}elsif((grep {$_=~ /\*/} @aa2) or (($reverse_start_codon2 ne "ATG") and ($reverse_start_codon2 ne "GTG")) or (($reverse_stop_codon2 ne "TAA") and ($reverse_stop_codon2 ne "TAG") and ($reverse_stop_codon2 ne "TGA"))){#non-standard start or stop codon for _CDS_aa
							my $gene_length=($start2-$end2+1);
							my ($start_left,$start_right,$end_left,$end_right);
							$start_left=$start2-60;
							if (($length-$start2)>=9000) {
								$start_right=$start2+9000;
							}elsif (($length-$start2)<9000){
								$start_right=$start2+(int(($length-$start2)/3)-1)*3;
							}
							if ($end2>=9000) {
								$end_left=$end2-9000;
							}elsif ($end2<9000){
								$end_left=$end2-(int($end2/3)-1)*3;
							}
							$end_right=$end2+int($gene_length/3)*3;

							my $str1=substr($sequence,($start_left),($start_right-$start_left));
							my $str2=substr($sequence,($end_left-1),($end_right-$end_left));
							my $length1=length $str1;
							my $length2=length $str2;
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $aa1;#left range for start codon
							#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com1,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa1.=$hash_codon{$codon};
								}else{
									$aa1.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa1=split //,$aa1;
							my (@star1,@start1);
							foreach (0..$#aa1){
								if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
									push @star1,$_;
								}
								if ($aa1[$_]=~ /M/){
									push @start1,$_;
								}
							}
							my $right_star_position=pop @star1;
							my @start2;
							foreach (@start1){
								push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
							}
							my $right_start_position=$start2[0];
							my $aa2;#right range for stop codon
							#for (my $i=0;$i<$length2;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length2-3);$i+=3){# delete stop codon
								my $codon=substr ($rev_com2,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa2.=$hash_codon{$codon};
								}else{
									$aa2.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa2=split //,$aa2;
							my @star2;
							foreach (0..$#aa2){
								if ($aa2[$_]=~ /\*/){
									push @star2,$_;
								}
							}
							my $left_star_position=shift @star2;
							if ((defined $right_start_position) and (defined $left_star_position)){
								my $start=$start_right-$right_start_position * 3;
								my $end=$end_right-($left_star_position * 3+3);
								my $str=substr($sequence,($end-1),($start-$end+1));
								my $length=length $str;
								my $rev_com=reverse $str;
								$rev_com=~ tr/ACGTacgt/TGCAtgca/;
								my $aa;
								#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa.=$hash_codon{$codon};
									}else{
										$aa.="X";
										my $j=$i+1;
										print "Bad codon $codon in position $j of gene $name in species $contig!\n";
									}
								}
								print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
							}elsif ((!defined $right_start_position) and (defined $left_star_position)){
								my $end=$end_right-($left_star_position * 3+3);
								print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(".$end."..".$start2.")\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								print $logfile "Warning: $name (negative no-intron PCG) has alternative start codon!\n";
							}
						}
					}
				}
			}
		}

		############################################################
		## one- or two-intron protein-coding-gene and tRNA gene
		############################################################
		if((keys %{$hash{$name}} >= 2) and (keys %{$hash{$name}} <= 5)){
			my (@position1,@position2,@position3,@position4);
			my $coding;
			foreach my $gene_coding_CDS (sort {$b cmp $a} keys %{$hash{$name}}){
				if ($gene_coding_CDS=~ /^((.+)_gene)$/){
					$coding=$2;
					foreach (keys %{$hash{$name}{$1}}){
						push @position1,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /^(((.+)-1)_coding)$/ or $gene_coding_CDS=~ /^(((.+)-1)_coding_aa)$/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position2,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /^(((.+)-2)_coding)$/ or $gene_coding_CDS=~ /^(((.+)-2)_coding_aa)$/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position3,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}elsif($gene_coding_CDS=~ /^(((.+)-3)_coding)$/ or $gene_coding_CDS=~ /^(((.+)-3)_coding_aa)$/){
					foreach (keys %{$hash{$name}{$1}}){
						push @position4,$_."\t".$hash{$name}{$gene_coding_CDS}{$_};
					}
				}
			}

			while (@position1 and (@position2 or @position3)){# _gene
				my $coding1=$coding."-1_coding";
				my $coding2=$coding."-2_coding";
				my $coding3=$coding."-3_coding";

				my $position1=shift @position1;
				my ($start1,$end1,$number1)=(split /\t/,$position1)[0,1,2];# _gene
				my $position2=shift @position2;
				my ($start2,$end2,$number2);# -1_coding or -1_coding_aa
				my $position3=shift @position3;
				my ($start3,$end3,$number3);# -2_coding or -2_coding_aa
				my $position4=shift @position4;
				my ($start4,$end4,$number4)=(split /\t/,$position4)[0,1,2] if (defined $position4 ne "");# -3_coding or -3_coding_aa

				my (@coding1,@coding2);
				if (defined $position2 ne ""){
	    			($start2,$end2,$number2)=(split /\t/,$position2)[0,1,2];
				}elsif(defined $position2 eq ""){# very short exon (-1_coding)
					foreach (keys %hash_exon){
						if ($_ eq $coding1){
							foreach (sort {$hash_exon{$_}{$b} <=> $hash_exon{$_}{$a}} keys %{$hash_exon{$_}}){
								push @coding1,$_;
							}
							$start2=$start1;
							if ($start1 < $end1){
								$end2=$start1+$coding1[0]-1;
							}elsif($start1 > $end1){
								$end2=$start1-$coding1[0]+1;
							}
						}
					}
				}
				if (defined $position3 ne ""){
					($start3,$end3,$number3)=(split /\t/,$position3)[0,1,2];
				}elsif(defined $position3 eq ""){# very short exon (-2_coding)
					foreach (keys %hash_exon){
						if ($_ eq $coding2){
							foreach (sort {$hash_exon{$_}{$b} <=> $hash_exon{$_}{$a}} keys %{$hash_exon{$_}}){
								push @coding2,$_;
							}
							$start3=$start1;
							if ($start1 < $end1){
								$end3=$start1+$coding2[0]-1;
							}elsif($start1 > $end1){
								$end3=$start1-$coding2[0]+1;
							}
						}
					}
				}

				############################################################
				## one-intron tRNA
				############################################################
				if($name=~ /(\D+)-(\D){3}/){
					if (($start1 == $start2) and ($end1 == $end3)){# identical tRNA boundary for _gene and -1_coding, -2_coding
						if ($start1 < $end1){
							print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."tRNA"."            "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
						}elsif($start1 > $end1){
							print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."tRNA"."            "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
						}
					}elsif(($start1 != $start2) or ($end1 != $end3)){# non-identical tRNA boundary for _gene and -1_coding, -2_coding
						if ($start1 < $end1){
							print $out_annotation "     "."gene"."            ".$start2."..".$end3."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."tRNA"."            "."join(".$start2."..".$end2.",".$start3."..".$end3.")"."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $logfile "Warning: $name (positive one-intron tRNA) has non-identical start or end boundary!\n";
						}elsif($start1 > $end1){
							print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start2.")"."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "     "."tRNA"."            "."complement(join(".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
							print $out_annotation "                     "."/gene=\"$name\""."\n";
							print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
							#print $logfile "Warning: $name (negative one-intron tRNA) has non-identical start or end boundary!\n";
						}
					}

				############################################################
				## one- or two-intron protein-coding-gene
				############################################################
				}elsif ($name!~ /(\D+)-(\D){3}/){
					#print "$name\n$start1\t$end1\t$start2\t$end2\t$start3\t$end3\t$start4\t$end4\n";
					my $length_ref_exon1=$hash_exon_length{"$name-1_coding"};
					my $length_ref_exon2=$hash_exon_length{"$name-2_coding"};
					my $sequence_ref_exon1=$hash_exon_sequence{"$name-1_coding"};
					my $sequence_ref_exon2=$hash_exon_sequence{"$name-2_coding"};

					############################################################
					## one-intron protein-coding-gene
					############################################################
					if (($start1 == $start2) and ($end1 == $end3)){# identical PCG boundary for _gene and -1_coding_aa, -2_coding_aa
						if ($start1 < $end1){# positive
							my $str=substr($sequence,($start1-1),($end1-$start1+1));
							my $length=length $str;
							my $forward_start_codon=substr($str,0,3);
							my $forward_stop_codon=substr($str,($length-3),3);
							my $str1=substr($sequence,($start1-1),($end2-$start1+1));
							my $str2=substr($sequence,($start3-1),($end1-$start3+1));
							my $str3=$str1.$str2;
							my $length_exon=length $str3;

							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left;
									if ($start1>=9000) {
										$start_left=$start1-9000;
									}elsif ($start1<9000){
										$start_left=$start1-(int($start1/3)-1)*3;
									}
									my $start_right=$start1+60;
									my $length_exon1=($end2-$start1+1);
									my $end_right;
									if (($length-$start3)>=9000) {
										$end_right=$start3+9000;
									}elsif (($length-$start3)<9000){
										$end_right=$start3+((int(($length-$start3)/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
									my $str2=substr($sequence,($start3-1),($end_right-$start3));
									my $length1=length $str1;
									my $length2=length $str2;
									my $aa1;#left range for start codon
									#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($str1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if ($aa1[$_]=~ /M/){
											push @start1,$_;
											}
									}
									my $left_star_position=pop @star1;
									until ($left_star_position<(3000+$length_exon1/3)){
										$left_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
									}
									my $left_start_position=$start2[0];
									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($str2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $right_star_position=shift @star2;
									if ((defined $left_start_position) and (defined $right_star_position)){
										my $start=$start_left+$left_start_position * 3;
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2+1;
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2+2;
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start-1),($end2_new-($start-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $left_start_position) and (defined $right_star_position)){
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2+1;
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2+2;
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start1."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (positive one-intron PCG) has alternative start codon!\n" if ($name ne "rps12+2");
									}else{
										print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (positive one-intron PCG) maybe pseudogene!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left;
									if ($start1>=9000) {
										$start_left=$start1-9000;
									}elsif ($start1<9000){
										$start_left=$start1-(int($start1/3)-1)*3;
									}
									my $start_right=$start1+60;
									my $length_exon1=($end2-$start1+1);
									my $end_right;
									if (($length-$start3)>=9000) {
										$end_right=$start3+9000;
									}elsif (($length-$start3)<9000){
										$end_right=$start3+((int(($length-$start3)/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
									my $str2=substr($sequence,($start3-1),($end_right-$start3));
									my $length1=length $str1;
									my $length2=length $str2;
									my $coding_1_length=$hash_exon_length{"$name-1_coding"};
									my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
									my $left_start_position;
									for (my $i=0;$i<($length1-3);$i+=1){
										my $exon=substr ($str1,$i,$coding_1_length);
										$exon=lc $exon;
										if ($exon eq $coding_1_sequence) {
											$left_start_position=$i;
										}
									}

									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($str2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $right_star_position=shift @star2;
									if ((defined $left_start_position) and (defined $right_star_position)){
										my $start=$start_left+$left_start_position;
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										$end2_new=$start+$coding_1_length-1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start-1),($end2_new-($start-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $left_start_position) and (defined $right_star_position)){
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										$end2_new=$start1+$coding_1_length-1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start1."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (positive one-intron PCG) has alternative start codon!\n";
									}
								}
							}
						}elsif($start1 > $end1){# negative
							my $str=substr($sequence,($end1-1),($start1-$end1+1));
							my $rev_com=reverse $str;
							$rev_com=~ tr/ACGTacgt/TGCAtgca/;
							my $length=length $str;
							my $reverse_start_codon=substr($rev_com,0,3);
							my $reverse_stop_codon=substr($rev_com,($length-3),3);
							my $str1=substr($sequence,($end2-1),($start1-$end2+1));
							my $str2=substr($sequence,($end1-1),($start3-$end1+1));
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $str3=$rev_com1.$rev_com2;
							my $length_exon=length $str3;
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left=$start1-60;
									my $start_right;
									if (($length-$start1)>=9000) {
										$start_right=$start1+9000;
									}elsif (($length-$start1)<9000){
										$start_right=$start1+(int(($length-$start1)/3)-1)*3;
									}
									my $length_exon1=($start1-$end2+1);
									my $end_left;
									if ($start3>=9000) {
										$end_left=$start3-9000;
									}elsif ($start3<9000){
										$end_left=$start3-((int($start3/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left),($start_right-$start_left));
									my $str2=substr($sequence,($end_left),($start3-$end_left));
									my $length1=length $str1;
									my $length2=length $str2;
									my $rev_com1=reverse $str1;
									$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_com2=reverse $str2;
									$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa1;#left range for start codon
									#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if ($aa1[$_]=~ /M/){
											push @start1,$_;
										}
									}
									my $right_star_position=pop @star1;
									until ($right_star_position<(3000+$length_exon1/3)){
										$right_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
									}
									my $right_start_position=$start2[0];
									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($rev_com2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $left_star_position=shift @star2;
									if ((defined $right_start_position) and (defined $left_star_position)){
										my $start=$start_right-$right_start_position * 3;
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2-1;
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2-2;
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $right_start_position) and (defined $left_star_position)){
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2-1;
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2-2;
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;
										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start1."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (negative one-intron PCG) has alternative start codon!\n" if ($name ne "rps12+2");
									}else{
										print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (negative one-intron PCG) maybe pseudogene!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left=$start1-60;
									my $start_right;
									if (($length-$start1)>=9000) {
										$start_right=$start1+9000;
									}elsif (($length-$start1)<9000){
										$start_right=$start1+(int(($length-$start1)/3)-1)*3;
									}
									my $length_exon1=($start1-$end2+1);
									my $end_left;
									if ($start3>=9000) {
										$end_left=$start3-9000;
									}elsif ($start3<9000){
										$end_left=$start3-((int($start3/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left),($start_right-$start_left));
									my $str2=substr($sequence,($end_left),($start3-$end_left));
									my $length1=length $str1;
									my $length2=length $str2;
									my $rev_com1=reverse $str1;
									$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_com2=reverse $str2;
									$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
									my $coding_1_length=$hash_exon_length{"$name-1_coding"};
									my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
									my $right_start_position;
									for (my $i=0;$i<($length1-3);$i+=1){
										my $exon=substr ($rev_com1,$i,$coding_1_length);
										$exon=lc $exon;
										if ($exon eq $coding_1_sequence) {
											$right_start_position=$i;
										}
									}

									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($rev_com2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $left_star_position=shift @star2;
									if ((defined $right_start_position) and (defined $left_star_position)){
										my $start=$start_right-$right_start_position;
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										$end2_new=$start-$coding_1_length+1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $right_start_position) and (defined $left_star_position)){
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										$end2_new=$start1-$coding_1_length+1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start1."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (negative one-intron PCG) has alternative start codon!\n";
									}
								}
							}
						}
					}elsif ((($start1 != $start2) and (!defined $start4)) or ((!defined $start4) and ($end1 != $end3))){# non-identical PCG boundary for _gene and -1_coding_aa, -2_coding_aa
						if ($start1 < $end1){# positive
							my $str1=substr($sequence,($start2-1),($end2-$start2+1));
							my $str2=substr($sequence,($start3-1),($end3-$start3+1));
							my $str3=$str1.$str2;
							my $length_exon=length $str3;
							my $forward_start_codon=substr($str3,0,3);
							my $forward_stop_codon=substr($str3,($length_exon-3),3);
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            ".$start2."..".$end3."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2.",".$start3."..".$end3.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left;
									if ($start2>=9000) {
										$start_left=$start2-9000;
									}elsif ($start2<9000){
										$start_left=$start2-(int($start2/3)-1)*3;
									}
									my $start_right=$start2+60;
									my $length_exon1=($end2-$start2+1);
									my $end_right;
									if (($length-$start3)>=9000) {
										$end_right=$start3+9000;
									}elsif (($length-$start3)<9000){
										$end_right=$start3+((int(($length-$start3)/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left-1),($start_right-$start_left));
									my $str2=substr($sequence,($start3-1),($end_right-$start3));
									my $length1=length $str1;
									my $length2=length $str2;
									my $aa1;#left range for start codon
									#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($str1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if ($aa1[$_]=~ /M/){
											push @start1,$_;
											}
									}
									my $left_star_position=pop @star1;
									until ($left_star_position<(3000+$length_exon1/3)){
										$left_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
									}
									my $left_start_position=$start2[0];
									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($str2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $right_star_position=shift @star2;
									if ((defined $left_start_position) and (defined $right_star_position)){
										my $start=$start_left+$left_start_position * 3;
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2+1;
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2+2;
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start-1),($end2_new-($start-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $left_start_position) and (defined $right_star_position)){
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2+1;
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2+2;
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (positive one-intron PCG) has alternative start codon!\n" if ($name ne "rps12+2");
									}else{
										print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end1.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (positive one-intron PCG) maybe pseudogene!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            ".$start2."..".$end3."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2.",".$start3."..".$end3.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left;
									if ($start2>=9000) {
										$start_left=$start2-9000;
									}elsif ($start2<9000){
										$start_left=$start2-(int($start2/3)-1)*3;
									}
									my $start_right=$start2+60;
									my $length_exon1=($end2-$start2+1);
									my $end_right;
									if (($length-$start3)>=9000) {
										$end_right=$start3+9000;
									}elsif (($length-$start3)<9000){
										$end_right=$start3+((int(($length-$start3)/3)-1)*3);
									}
									my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
									my $str2=substr($sequence,($start3-1),($end_right-$start3));
									my $length1=length $str1;
									my $length2=length $str2;
									my $coding_1_length=$hash_exon_length{"$name-1_coding"};
									my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
									my $left_start_position;
									for (my $i=0;$i<($length1-3);$i+=1){
										my $exon=substr ($str1,$i,$coding_1_length);
										$exon=lc $exon;
										if ($exon eq $coding_1_sequence) {
											$left_start_position=$i;
										}
									}

									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($str2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $right_star_position=shift @star2;
									if ((defined $left_start_position) and (defined $right_star_position)){
										my $start=$start_left+$left_start_position;
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										$end2_new=$start+$coding_1_length-1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start-1),($end2_new-($start-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $left_start_position) and (defined $right_star_position)){
										my $end=($start3-1)+($right_star_position * 3+3);

										my ($str1,$str2,$str,$length,$end2_new,$start3_new);
										$end2_new=$start1+$coding_1_length-1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3-2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3-1;
										}
										$str1=substr($sequence,($start1-1),($end2_new-($start1-1)));
										$str2=substr($sequence,($start3_new-1),($end-($start3_new-1)));
										$str=$str1.$str2;
										$length=length $str;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($str,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2_new.",".$start3_new."..".$end.")"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (positive one-intron PCG) has alternative start codon!\n";
									}
								}
							}
						}elsif($start1 > $end1){# negative
							my $str1=substr($sequence,($end2-1),($start2-$end2+1));
							my $str2=substr($sequence,($end3-1),($start3-$end3+1));
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $str3=$rev_com1.$rev_com2;
							my $length_exon=length $str3;
							my $reverse_start_codon=substr($str3,0,3);
							my $reverse_stop_codon=substr($str3,($length_exon-3),3);
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str3,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($name ne "rpl16") and ($name ne "petB") and ($name ne "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start2.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left=$start2-60;
									my $start_right;
									if (($length-$start2)>=9000) {
										$start_right=$start2+9000;
									}elsif (($length-$start2)<9000){
										$start_right=$start2+(int(($length-$start2)/3)-1)*3;
									}
									my $length_exon1=($start2-$end2+1);
									my $end_left;
									if ($start3>=9000) {
										$end_left=$start3-9000;
									}elsif ($start3<9000){
										$end_left=$start3-(int($start3/3)-1)*3;
									}
									my $str1=substr($sequence,($start_left),($start_right-$start_left));
									my $str2=substr($sequence,($end_left),($start3-$end_left));
									my $length1=length $str1;
									my $length2=length $str2;
									my $rev_com1=reverse $str1;
									$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_com2=reverse $str2;
									$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
									my $aa1;#left range for start codon
									#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com1,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa1.=$hash_codon{$codon};
										}else{
											$aa1.="X";
											my $m=$i+1;
											print "Bad codon $codon in position $m of gene $name in species $contig!\n";
										}
									}
									my @aa1=split //,$aa1;
									my (@star1,@start1);
									foreach (0..$#aa1){
										if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
											push @star1,$_;
										}
										if ($aa1[$_]=~ /M/){
											push @start1,$_;
										}
									}
									my $right_star_position=pop @star1;
									until ($right_star_position<(3000+$length_exon1/3)){
										$right_star_position=pop @star1;
									}
									my @start2;
									foreach (@start1){
										push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
									}
									my $right_start_position=$start2[0];
									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($rev_com2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $left_star_position=shift @star2;
									if ((defined $right_start_position) and (defined $left_star_position)){
										my $start=$start_right-$right_start_position * 3;
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2-1;
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2-2;
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $right_start_position) and (defined $left_star_position)){
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										if ($length_ref_exon1 % 3==0) {
											$end2_new=$end2;
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$end2_new=$end2-1;
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$end2_new=$end2-2;
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start2."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (negative one-intron PCG) has alternative start codon!\n" if ($name ne "rps12+2");
									}else{
										print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start3.",".$end2."..".$start1."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (negative one-intron PCG) maybe pseudogene!\n";
									}
								}
							}elsif (($name eq "rpl16") or ($name eq "petB") or ($name eq "petD")) {
								if (($length_ref_exon1 % 3==0) and ($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon for _gene
									print $out_annotation "     "."gene"."            "."complement(".$end3."..".$start2.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
								}elsif(($length_ref_exon1 % 3!=0) or (grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
									my $start_left=$start2-60;
									my $start_right;
									if (($length-$start2)>=9000) {
										$start_right=$start2+9000;
									}elsif (($length-$start2)<9000){
										$start_right=$start2+(int(($length-$start2)/3)-1)*3;
									}
									my $length_exon1=($start2-$end2+1);
									my $end_left;
									if ($start3>=9000) {
										$end_left=$start3-9000;
									}elsif ($start3<9000){
										$end_left=$start3-(int($start3/3)-1)*3;
									}
									my $str1=substr($sequence,($start_left),($start_right-$start_left));
									my $str2=substr($sequence,($end_left),($start3-$end_left));
									my $length1=length $str1;
									my $length2=length $str2;
									my $rev_com1=reverse $str1;
									$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
									my $rev_com2=reverse $str2;
									$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
									my $coding_1_length=$hash_exon_length{"$name-1_coding"};
									my $coding_1_sequence=$hash_exon_sequence{"$name-1_coding"};
									my $right_start_position;
									for (my $i=0;$i<($length1-3);$i+=1){
										my $exon=substr ($rev_com1,$i,$coding_1_length);
										$exon=lc $exon;
										if ($exon eq $coding_1_sequence) {
											$right_start_position=$i;
										}
									}

									my $aa2;#right range for stop codon
									#for (my $k=0;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $k=0;$k<($length2-3);$k+=3){# delete stop codon
										my $codon=substr ($rev_com2,$k,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa2.=$hash_codon{$codon};
										}else{
											$aa2.="X";
											my $n=$k+1;
											print "Bad codon $codon in position $n of gene $name in species $contig!\n";
										}
									}
									my @aa2=split //,$aa2;
									my @star2;
									foreach (0..$#aa2){
										if ($aa2[$_]=~ /\*/){
											push @star2,$_;
										}
									}
									my $left_star_position=shift @star2;
									if ((defined $right_start_position) and (defined $left_star_position)){
										my $start=$start_right-$right_start_position;
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										$end2_new=$start-$coding_1_length+1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									}elsif ((!defined $right_start_position) and (defined $left_star_position)){
										my $end=($start3+1)-($left_star_position * 3+3);

										my ($str1,$str2,$str,$length,$rev_com,$end2_new,$start3_new);
										$end2_new=$start1-$coding_1_length+1;
										if ($length_ref_exon1 % 3==0) {
											$start3_new=$start3;
										}elsif ($length_ref_exon1 % 3==1) {
											$start3_new=$start3+2;
										}elsif ($length_ref_exon1 % 3==2) {
											$start3_new=$start3+1;
										}
										$str1=substr($sequence,($end2_new-1),($start1-($end2_new-1)));
										$str2=substr($sequence,($end-1),($start3_new-($end-1)));
										$str=$str1.$str2;
										$length=length $str;
										$rev_com=reverse $str;
										$rev_com=~ tr/ACGTacgt/TGCAtgca/;

										my $aa;
										#for (my $i=0;$i<$length;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
										for (my $i=0;$i<($length-3);$i+=3){# delete stop codon
											my $codon=substr ($rev_com,$i,3);
											$codon=uc $codon;
											if (exists $hash_codon{$codon}){
												$aa.=$hash_codon{$codon};
											}else{
												$aa.="X";
												my $j=$i+1;
												print "Bad codon $codon in position $j of gene $name in species $contig!\n";
											}
										}
										print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start3_new.",".$end2_new."..".$start2."))"."\n";
										print $out_annotation "                     "."/gene=\"$name\""."\n";
										print $out_annotation "                     "."/codon_start=1"."\n";
										print $out_annotation "                     "."/transl_table=11"."\n";
										print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
										#print $out_annotation "                     "."/translation=\"$aa\""."\n";
										print $logfile "Warning: $name (negative one-intron PCG) has alternative start codon!\n";
									}
								}
							}
						}
					}

					############################################################
					## two-intron protein-coding-gene
					############################################################
					if(($start1 == $start2) and (defined $end4) and ($end1 == $end4)){# identical PCG boundary for _gene and -1_coding, -3_coding
						if ($start1 < $end1){# positive
							my $str=substr($sequence,($start1-1),($end1-$start1+1));
							my $length_3=length $str;
							my $forward_start_codon=substr($str,0,3);
							my $forward_stop_codon=substr($str,($length_3-3),3);
							my $str1=substr($sequence,($start1-1),($end2-$start1+1));
							my $str2=substr($sequence,($start3-1),($end3-$start3+1));
							my $str3=substr($sequence,($start4-1),($end1-$start4+1));
							my $str4=$str1.$str2.$str3;
							my $length_exon=length $str4;
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
								print $out_annotation "     "."gene"."            ".$start1."..".$end1."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end3.",".$start4."..".$end1.")"."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left;
								if ($start1>=9000) {
									$start_left=$start1-9000;
								}elsif ($start1<9000){
									$start_left=$start1-int($start1/3)*3;
								}
								my $start_right=$start1+60;
								my $length_exon1=($end2-$start1+1);
								my $length_exon2=($end3-$start3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_right;
								if (($length-$start4)>=9002) {
									if ($length_exon12 % 3==0){
										$end_right=$start4+9000;
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+9001;
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+9002;
									}
								}elsif (($length-$start4)<9002){
									if ($length_exon12 % 3==0){
										$end_right=$start4+(int(($length-$start4)/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+(int(($length-$start4)/3)*3+1);
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+(int(($length-$start4)/3)*3+2);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($start4-1),($end_right-$start4+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $aa1;#left range for start codon
								#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($str1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if ($aa1[$_]=~ /M/){
										push @start1,$_;
										}
								}
								my $left_star_position=pop @star1;
								until ($left_star_position<(3000+$length_exon1/3)){
									$left_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
								}
								my $left_start_position=$start2[0];
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								#for (my $k=$j;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($str2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $right_star_position=shift @star2;
								if ((defined $left_start_position) and (defined $right_star_position)){
									my $start=$start_left+$left_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($left_star_position * 3+3);
									}
									my $str1=substr($sequence,($start-1),($end2-$start+1));
									my $str2=substr($sequence,($start3-1),($end3-$start3+1));
									my $str3=substr($sequence,($start4-1),($end-$start4+1));
									my $str=$str1.$str2.$str3;
									my $length_exon=length $str;
									my $aa;
									#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
										my $codon=substr ($str,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2.",".$start3."..".$end3.",".$start4."..".$end.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								}elsif ((!defined $left_start_position) and (defined $right_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($left_star_position * 3+3);
									}
									print $out_annotation "     "."gene"."            ".$start1."..".$end."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start1."..".$end2.",".$start3."..".$end3.",".$start4."..".$end.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									print $logfile "Warning: $name (positive two-intron PCG) has alternative start codon!\n";
								}
							}else{
								print "The exon length of $name (positive two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}elsif($start1 > $end1){# negative
							my $str=substr($sequence,($end1-1),($start1-$end1+1));
							my $rev_com=reverse $str;
							$rev_com=~ tr/ACGTacgt/TGCAtgca/;
							my $length_3=length $str;
							my $reverse_start_codon=substr($rev_com,0,3);
							my $reverse_stop_codon=substr($rev_com,($length_3-3),3);
							my $str1=substr($sequence,($end2-1),($start1-$end2+1));
							my $str2=substr($sequence,($end3-1),($start3-$end3+1));
							my $str3=substr($sequence,($end1-1),($start4-$end1+1));
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com3=reverse $str3;
							$rev_com3=~ tr/ACGTacgt/TGCAtgca/;
							my $str4=$rev_com1.$rev_com2.$rev_com3;
							my $length_exon=length $str4;
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon
								print $out_annotation "     "."gene"."            "."complement(".$end1."..".$start1.")"."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(join(".$end1."..".$start4.",".$end3."..".$start3.",".$end2."..".$start1."))"."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left=$start1-60;
								my $start_right;
								if (($length-$start1)>=9000) {
									$start_right=$start1+9000;
								}elsif (($length-$start1)<9000){
									$start_right=$start1+int(($length-$start1)/3)*3;
								}
								my $length_exon1=($start1-$end2+1);
								my $length_exon2=($start3-$end3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_left;
								if ($start4>=9002) {
									if ($length_exon12 % 3==0){
										$end_left=$start4-9000;
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-9002;
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-9001;
									}
								}elsif ($start4<9002){
									if ($length_exon12 % 3==0){
										$end_left=$start4-(int($start4/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-(int($start4/3)*3+2);
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-(int($start4/3)*3+1);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($end_left-1),($start4-$end_left+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $rev_com1=reverse $str1;
								$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_com2=reverse $str2;
								$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
								my $aa1;#left range for start codon
								#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if ($aa1[$_]=~ /M/){
										push @start1,$_;
									}
								}
								my $right_star_position=pop @star1;
								until ($right_star_position<(3000+$length_exon1/3)){
									$right_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
								}
								my $right_start_position=$start2[0];
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								#for (my $k=$j;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($rev_com2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $left_star_position=shift @star2;
								if ((defined $right_start_position) and (defined $left_star_position)){
									my $start=$start_right-$right_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}
									my $str1=substr($sequence,($end2-1),($start-$end2+1));
									my $str2=substr($sequence,($end3-1),($start3-$end3+1));
									my $str3=substr($sequence,($end-1),($start4-$end+1));
									my $str=$str1.$str2.$str3;
									my $length_exon=length $str;
									my $rev_com=reverse $str;
									$rev_com=~ tr/ACGTacgt/TGCAtgca/;
									my $aa;
									#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4.",".$end3."..".$start3.",".$end2."..".$start."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								}elsif ((!defined $right_start_position) and (defined $left_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}
									print $out_annotation "     "."gene"."            "."complement(".$end."..".$start1.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4.",".$end3."..".$start3.",".$end2."..".$start1."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									print $logfile "Warning: $name (negative two-intron PCG) has alternative start codon!\n";
								}
							}else{
								print "The exon length of $name (negative two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}
					}elsif((($start1 != $start2) and (defined $end4)) or ((defined $end4) and ($end1 != $end4))){# non-identical PCG boundary for _gene and -1_coding, -3_coding
						if ($start1 < $end1){# positive
							my $str1=substr($sequence,($start2-1),($end2-$start2+1));
							my $str2=substr($sequence,($start3-1),($end3-$start3+1));
							my $str3=substr($sequence,($start4-1),($end4-$start4+1));
							my $str4=$str1.$str2.$str3;
							my $length_exon=length $str4;
							my $forward_start_codon=substr($str4,0,3);
							my $forward_stop_codon=substr($str4,($length_exon-3),3);
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($forward_start_codon eq "ATG") or ($forward_start_codon eq "GTG")) and (($forward_stop_codon eq "TAA") or ($forward_stop_codon eq "TAG") or ($forward_stop_codon eq "TGA"))){# standard start and stop codon for _gene
								print $out_annotation "     "."gene"."            ".$start2."..".$end4."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2.",".$start3."..".$end3.",".$start4."..".$end4.")"."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($forward_start_codon ne "ATG") and ($forward_start_codon ne "GTG")) or (($forward_stop_codon ne "TAA") and ($forward_stop_codon ne "TAG") and ($forward_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left;
								if ($start2>=9000) {
									$start_left=$start2-9000;
								}elsif ($start2<9000){
									$start_left=$start2-int($start2/3)*3;
								}
								my $start_right=$start2+60;
								my $length_exon1=($end2-$start2+1);
								my $length_exon2=($end3-$start3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_right;
								if (($length-$start4)>=9002) {
									if ($length_exon12 % 3==0){
										$end_right=$start4+9000;
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+9001;
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+9002;
									}
								}elsif (($length-$start4)<9002){
									if ($length_exon12 % 3==0){
										$end_right=$start4+(int(($length-$start4)/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_right=$start4+(int(($length-$start4)/3)*3+1);
									}
									if ($length_exon12 % 3==2){
										$end_right=$start4+(int(($length-$start4)/3)*3+2);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($start4-1),($end_right-$start4+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $aa1;#left range for start codon
								#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($str1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if ($aa1[$_]=~ /M/){
										push @start1,$_;
										}
								}
								my $left_star_position=pop @star1;
								until ($left_star_position<(3000+$length_exon1/3)){
									$left_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $left_star_position) and ($_ > $left_star_position))
								}
								my $left_start_position=$start2[0];
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								#for (my $k=$j;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($str2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $right_star_position=shift @star2;
								if ((defined $left_start_position) and (defined $right_star_position)){
									my $start=$start_left+$left_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($left_star_position * 3+3);
									}
									my $str1=substr($sequence,($start-1),($end2-$start+1));
									my $str2=substr($sequence,($start3-1),($end3-$start3+1));
									my $str3=substr($sequence,($start4-1),($end-$start4+1));
									my $str=$str1.$str2.$str3;
									my $length_exon=length $str;
									my $aa;
									#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
										my $codon=substr ($str,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									print $out_annotation "     "."gene"."            ".$start."..".$end."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start."..".$end2.",".$start3."..".$end3.",".$start4."..".$end.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								}elsif ((!defined $left_start_position) and (defined $right_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4-1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4+1)+($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)+($left_star_position * 3+3);
									}
									print $out_annotation "     "."gene"."            ".$start2."..".$end."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."join(".$start2."..".$end2.",".$start3."..".$end3.",".$start4."..".$end.")"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									print $logfile "Warning: $name (positive two-intron PCG) has alternative start codon!\n";
								}
							}else{
								print "The exon length of $name (positive two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}elsif($start1 > $end1){# negative
							my $str=substr($sequence,($end4-1),($start2-$end4+1));
							my $rev_com=reverse $str;
							$rev_com=~ tr/ACGTacgt/TGCAtgca/;
							my $length_3=length $str;
							my $reverse_start_codon=substr($rev_com,0,3);
							my $reverse_stop_codon=substr($rev_com,($length_3-3),3);
							my $str1=substr($sequence,($end2-1),($start2-$end2+1));
							my $str2=substr($sequence,($end3-1),($start3-$end3+1));
							my $str3=substr($sequence,($end4-1),($start4-$end4+1));
							my $rev_com1=reverse $str1;
							$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com2=reverse $str2;
							$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
							my $rev_com3=reverse $str3;
							$rev_com3=~ tr/ACGTacgt/TGCAtgca/;
							my $str4=$rev_com1.$rev_com2.$rev_com3;
							my $length_exon=length $str4;
							my $aa_exon;
							#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
							for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
								my $codon=substr ($str4,$i,3);
								$codon=uc $codon;
								if (exists $hash_codon{$codon}){
									$aa_exon.=$hash_codon{$codon};
								}else{
									$aa_exon.="X";
									my $j=$i+1;
									print "Bad codon $codon in position $j of gene $name in species $contig!\n";
								}
							}
							my @aa_exon=split //,$aa_exon;
							if (($length_exon % 3==0) and (!grep {$_=~ /\*/} @aa_exon) and (($reverse_start_codon eq "ATG") or ($reverse_start_codon eq "GTG")) and (($reverse_stop_codon eq "TAA") or ($reverse_stop_codon eq "TAG") or ($reverse_stop_codon eq "TGA"))){# standard start and stop codon
								print $out_annotation "     "."gene"."            "."complement(".$end4."..".$start2.")"."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "     "."CDS"."             "."complement(join(".$end4."..".$start4.",".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
								print $out_annotation "                     "."/gene=\"$name\""."\n";
								print $out_annotation "                     "."/codon_start=1"."\n";
								print $out_annotation "                     "."/transl_table=11"."\n";
								print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
								#print $out_annotation "                     "."/translation=\"$aa_exon\""."\n";
							}elsif((grep {$_=~ /\*/} @aa_exon) or (($reverse_start_codon ne "ATG") and ($reverse_start_codon ne "GTG")) or (($reverse_stop_codon ne "TAA") and ($reverse_stop_codon ne "TAG") and ($reverse_stop_codon ne "TGA"))){# non-standard start or stop codon for _gene
								my $start_left=$start2-60;
								my $start_right;
								if (($length-$start2)>=9000) {
									$start_right=$start2+9000;
								}elsif (($length-$start2)<9000){
									$start_right=$start2+int(($length-$start2)/3)*3;
								}
								my $length_exon1=($start2-$end2+1);
								my $length_exon2=($start3-$end3+1);
								my $length_exon12=$length_exon1+$length_exon2;
								my $end_left;
								if ($start4>=9002) {
									if ($length_exon12 % 3==0){
										$end_left=$start4-9000;
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-9002;
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-9001;
									}
								}elsif ($start4<9002){
									if ($length_exon12 % 3==0){
										$end_left=$start4-(int($start4/3)*3);
									}
									if ($length_exon12 % 3==1){
										$end_left=$start4-(int($start4/3)*3+2);
									}
									if ($length_exon12 % 3==2){
										$end_left=$start4-(int($start4/3)*3+1);
									}
								}
								my $str1=substr($sequence,($start_left-1),($start_right-$start_left+1));
								my $str2=substr($sequence,($end_left-1),($start4-$end_left+1));
								my $length1=length $str1;
								my $length2=length $str2;
								my $rev_com1=reverse $str1;
								$rev_com1=~ tr/ACGTacgt/TGCAtgca/;
								my $rev_com2=reverse $str2;
								$rev_com2=~ tr/ACGTacgt/TGCAtgca/;
								my $aa1;#left range for start codon
								#for (my $i=0;$i<$length1;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $i=0;$i<($length1-3);$i+=3){# delete stop codon
									my $codon=substr ($rev_com1,$i,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa1.=$hash_codon{$codon};
									}else{
										$aa1.="X";
										my $m=$i+1;
										print "Bad codon $codon in position $m of gene $name in species $contig!\n";
									}
								}
								my @aa1=split //,$aa1;
								my (@star1,@start1);
								foreach (0..$#aa1){
									if (($aa1[$_]=~ /\*/) and ($_ < 3000)){
										push @star1,$_;
									}
									if ($aa1[$_]=~ /M/){
										push @start1,$_;
									}
								}
								my $right_star_position=pop @star1;
								until ($right_star_position<(3000+$length_exon1/3)){
									$right_star_position=pop @star1;
								}
								my @start2;
								foreach (@start1){
									push @start2,$_ if ((defined $right_star_position) and ($_ > $right_star_position))
								}
								my $right_start_position=$start2[0];
								my $aa2;#right range for stop codon
								my $j;
								if ($length_exon1 % 3==0){
									$j=0;
								}
								if ($length_exon1 % 3==1){
									$j=2;
								}
								if ($length_exon1 % 3==2){
									$j=1;
								}
								#for (my $k=$j;$k<$length2;$k+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
								for (my $k=$j;$k<($length2-3);$k+=3){# delete stop codon
									my $codon=substr ($rev_com2,$k,3);
									$codon=uc $codon;
									if (exists $hash_codon{$codon}){
										$aa2.=$hash_codon{$codon};
									}else{
										$aa2.="X";
										my $n=$k+1;
										print "Bad codon $codon in position $n of gene $name in species $contig!\n";
									}
								}
								my @aa2=split //,$aa2;
								my @star2;
								foreach (0..$#aa2){
									if ($aa2[$_]=~ /\*/){
										push @star2,$_;
									}
								}
								my $left_star_position=shift @star2;
								if ((defined $right_start_position) and (defined $left_star_position)){
									my $start=$start_right-$right_start_position * 3;
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}
									my $str1=substr($sequence,($end2-1),($start-$end2+1));
									my $str2=substr($sequence,($end3-1),($start3-$end3+1));
									my $str3=substr($sequence,($end-1),($start4-$end+1));
									my $str=$str1.$str2.$str3;
									my $length_exon=length $str;
									my $rev_com=reverse $str;
									$rev_com=~ tr/ACGTacgt/TGCAtgca/;
									my $aa;
									#for (my $i=0;$i<$length_exon;$i+=3){# remain stop codon, either (length ($seq)-1) or (length ($seq)-2) is OK
									for (my $i=0;$i<($length_exon-3);$i+=3){# delete stop codon
										my $codon=substr ($rev_com,$i,3);
										$codon=uc $codon;
										if (exists $hash_codon{$codon}){
											$aa.=$hash_codon{$codon};
										}else{
											$aa.="X";
											my $j=$i+1;
											print "Bad codon $codon in position $j of gene $name in species $contig!\n";
										}
									}
									print $out_annotation "     "."gene"."            "."complement(".$end."..".$start.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4.",".$end3."..".$start3.",".$end2."..".$start."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
								}elsif ((!defined $right_start_position) and (defined $left_star_position)){
									my $end;
									if ($length_exon12 % 3==0){
										$end=($start4+1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==1){
										$end=($start4-1)-($left_star_position * 3+3);
									}
									if ($length_exon12 % 3==2){
										$end=($start4)-($left_star_position * 3+3);
									}
									print $out_annotation "     "."gene"."            "."complement(".$end."..".$start2.")\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "     "."CDS"."             "."complement(join(".$end."..".$start4.",".$end3."..".$start3.",".$end2."..".$start2."))"."\n";
									print $out_annotation "                     "."/gene=\"$name\""."\n";
									print $out_annotation "                     "."/codon_start=1"."\n";
									print $out_annotation "                     "."/transl_table=11"."\n";
									print $out_annotation "                     "."/product=\"".$hash_product{$name}."\""."\n";
									#print $out_annotation "                     "."/translation=\"$aa\""."\n";
									print $logfile "Warning: $name (negative two-intron PCG) need to be manually checked! Alternative start codon or pseudugene!\n";
								}
							}else{
								print "The exon length of $name (negative two-intron PCG) is not an intergral multiple of 3!\n";
							}
						}
					}
				}
			}
		}
	}
	print $out_annotation "ORIGIN      \n";
	print $out_annotation "$sequence\n";
	print $out_annotation "//\n";
	print $logfile "\n";
	close $out_annotation;
	close $logfile;
	my $now4=&gettime;
	print "\n$now4 || Finish annotating the $j sequence: $header.fasta!\n";
	%hash_codon=();
	%hash_product=();
}
%hash_in_all2=();
%hash_exon=();
%hash_exon_length=();
%hash_exon_sequence=();
unlink("reference1.fasta");
unlink("reference2.fasta");
unlink("reference3.fasta");
unlink("reference4.fasta");
my $now5=&gettime;
print "$now5 || Finish annotating all sequences and total elapsed time is: ",time-$now1," seconds!\n";


########################################
##subroutines
########################################
sub getdate {
	my ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday)=(localtime(time));
	my %hash_month=("01"=>"JAN","02"=>"FEB","03"=>"MAR","04"=>"APR","05"=>"MAY","06"=>"JUN","07"=>"JUL","08"=>"AUG","09"=>"SEP","10"=>"OCT","11"=>"NOV","12"=>"DEC");
	$year+=1900;
	$mon=($mon<9)?"0".($mon+1):($mon+1);
	$day=($day<10)?"0$day":$day;
	$hour=($hour<10)?"0$hour":$hour;
	$min=($min<10)?"0$min":$min;
	$sec=($sec<10)?"0$sec":$sec;

	my $now="$day-$hash_month{$mon}-$year";
}

sub gettime {
	my ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday)=(localtime(time));
	$year+=1900;
	$mon=($mon<9)?"0".($mon+1):($mon+1);
	$day=($day<10)?"0$day":$day;
	$hour=($hour<10)?"0$hour":$hour;
	$min=($min<10)?"0$min":$min;
	$sec=($sec<10)?"0$sec":$sec;

	my $now="$year.$mon.$day $hour:$min:$sec";
}

sub argument{
	my @options=("help|h","ref|r:s","seq|s:s","out|o:s","type|t:s","log|l:s");
	my %options;
	GetOptions(\%options,@options);
	exec ("pod2usage $0") if ((keys %options)==0 or $options{'h'} or $options{'help'});
	if(!exists $options{'ref'}){
		print "***ERROR: No reference directory is assigned!!!\n";
		exec ("pod2usage $0");
	}elsif(!exists $options{'seq'}){
		print "***ERROR: No sequence directory is assigned!!!\n";
		exec ("pod2usage $0");
	}
	return \%options;
}

sub default{
	my ($default_value,$option)=@_;
	if(exists $global_options->{$option}){
		return $global_options->{$option};
	}
	return $default_value;
}


__DATA__

=head1 NAME

    PGA.pl Plastid Genome Annotation

=head1 COPYRIGHT

    copyright (C) 2017 Xiao-Jian Qu

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    Plastid Genome Annotation

=head1 SYNOPSIS

    PGA.pl -r -s [-o -t -l]
    Copyright (C) 2017 Xiao-Jian Qu
    Please contact <quxiaojian@mail.kib.ac.cn>, if you have any bugs or questions.

    [-h -help]       help information.
    [-r -ref]        required: input directory name containing gb reference file(s) that from the same or close family,order etc.(default: reference)
    [-s -seq]        required: input directory name containing fasta sequence files(s) that you want to annotate.(default: sequence)
    [-o -out]        optional: output directory name containing annotated genebank(gb) file(s).(default: gb)
    [-t -type]       optional: circular or linear for fasta sequence files(s).(default: circular)
    [-l -log]        optional: log file name containing warning information for annotated genebank(gb) file(s).(default: warning)

=cut
