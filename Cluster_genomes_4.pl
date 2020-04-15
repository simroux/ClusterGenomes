#!/usr/bin/env perl
use strict;
use autodie;
use Getopt::Long;
my $h='';
my $fasta_file='';
my $th_id=95;
my $th_cover=85;
my $cmd='';
my $out='';
my $tag_fna=0;
my $n_threads=1;
my $nucmer_dir='';
GetOptions ('help' => \$h, 'h' => \$h, 'f=s'=>\$fasta_file, 'i=s'=>\$th_id, 'c=s'=>\$th_cover, 'd=s'=>\$nucmer_dir, 'nofna'=>\$tag_fna, 't=s'=>\$n_threads);
if ($h==1 || $fasta_file eq ""){ # If asked for help or did not set up any argument
	print "# Script to cluster new contigs
# Arguments :
-f : input fasta file
-d : path to the bin directory of MUMMER (e.g. tools/Mummer/4.0.0b2/bin/)
-c : % of the sequence covered (default 85)
-i : % of identity (default 95)
-nofna : do not generate a fasta file of the seeds
-t : number of threads (default: 1)
e.g. -c 80% -i 95% -> look for sequences similar at 95% ANI covering 80% or more of the shortest sequence\n";
	die "\n";
}
 if (!($nucmer_dir=~/\/$/)){$nucmer_dir.="/";}

$fasta_file=~/(.*)\.f.+/;
my $root=$1;

my $mummer_file=$root."-nucmer.out";
my $mummer_file_out=$root."-nucmer.out.coords";
my $out_file=$root."_".$th_id."-".$th_cover.".clstr";
my $out_file_fasta=$root."_".$th_id."-".$th_cover.".fna";
if (-e $mummer_file_out){
	print " mummer already run, no need to do it again\n";
}
else{
	print "comparing the contigs with mummer -> $mummer_file_out\n";
	$cmd="$nucmer_dir/nucmer --maxmatch --nooptimize $fasta_file $fasta_file -t $n_threads -p $mummer_file";
	print "Running nucmer : $cmd\n";
	$out=`$cmd`;
	print "$out\n";
	my $delta_file=$root."-nucmer.out.delta";
	$cmd="$nucmer_dir/show-coords $delta_file > $mummer_file_out";
	print "Running coords : $cmd\n";
	$out=`$cmd`;
	print "$out\n";
}


my %seq_len;
open my $fa,"<",$fasta_file;
my $c_c='';
while(<$fa>){
	chomp($_);
	if ($_=~/^>(\S+)/){$c_c=$1;}
	else{$seq_len{$c_c}+=length($_);}
}
close $fa;


my $cover_result=$root."-cover.csv";
my %store;
my %store_tmp;
my %check;
my %seen_code;
if (!(-e $cover_result)){
	print "parsing $mummer_file_out -> $cover_result\n";
	my $mummer_file_out_sorted=$mummer_file_out.".sorted";
	if (!(-e $mummer_file_out_sorted)){&sort($mummer_file_out,$mummer_file_out_sorted);}
	open my $s1,">",$cover_result;
	print $s1 "### Subject,Query,# Matches (estimated),Hit length,Subject_length,Query_length\n";
	open my $tsv,"<",$mummer_file_out_sorted;
	my $tag=0;
	my $last;
	while(<$tsv>){
		chomp($_);
		if ($_=~/\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=/){$tag=1}
		elsif($tag==1){
			my @tab=split(" ",$_);
			shift(@tab); # NOW HAVE TO REMOVE THE FIRST ELEMENT WE USE FOR SORTING
			if ($tab[11] eq $tab[12]){next;} # same genome
			if ($tab[12] eq ""){next;} # Not a relevant line
			if (!defined($seq_len{$tab[11]}) || !defined($seq_len{$tab[12]})){next;} # we don't have a length, the coords is probably taken from another step when we had more sequences, we skip
			if ($tab[12] ne $last){
				# calculate and store result from last query
				print "### $last\n";
				foreach my $hit (keys %store_tmp){
					print "Hit $hit vs last $last \n";
					my $real_cover=scalar(keys %{$store_tmp{$hit}});
					$store{$hit}{$last}{"match"}-=$store{$hit}{$last}{"cover"}-$real_cover;
					$store{$hit}{$last}{"cover"}=$real_cover;
					print $s1 $hit.",".$last.",".$store{$hit}{$last}{"match"}.",".$store{$hit}{$last}{"cover"}.",".$seq_len{$hit}.",".$seq_len{$last}."\n";
				}
				# re-initialize variables
				$last=$tab[12];
				%store_tmp=();
				if (defined($seen_code{$tab[12]})){
					print "!!! We ALREADY SAW $tab[12] ?!?!?!\n";
					<STDIN>;
				}
				$seen_code{$last}=1;
			}
			if ($tab[1]>$tab[0]){}
			else{
				my $temp=$tab[1];
				$tab[1]=$tab[0];
				$tab[0]=$temp;
			}
			if ($tab[4]>$tab[3]){}
			else{
				my $temp=$tab[4];
				$tab[4]=$tab[3];
				$tab[3]=$temp;
			}
			my $pcent=$tab[9];
			my $seed=$tab[11];
			my $seed_c=$tab[6];
			my $query=$tab[12];
			my $query_c=$tab[7];
			for(my $i=$tab[0];$i<=$tab[1];$i++){
				$store_tmp{$seed}{$i}++;
			}
			$store{$seed}{$query}{"match"}+=$seed_c*$pcent/100;
			$store{$seed}{$query}{"cover"}+=$seed_c;
		}
	}
	close $tsv;
	# We process the last one
	foreach my $hit (keys %store_tmp){
		my $real_cover=scalar(keys %{$store_tmp{$hit}});
		$store{$hit}{$last}{"match"}-=$store{$hit}{$last}{"cover"}-$real_cover;
		$store{$hit}{$last}{"cover"}=$real_cover;
		print $s1 $hit.",".$last.",".$store{$hit}{$last}{"match"}.",".$store{$hit}{$last}{"cover"}.",".$seq_len{$hit}.",".$seq_len{$last}."\n";
	}
	close $s1;
}
else{
	print "Reading $cover_result\n";
	open my $csv,"<",$cover_result;
	while(<$csv>){
		chomp($_);
		my @tab=split(",",$_);
		$store{$tab[0]}{$tab[1]}{"match"}=$tab[2];
		$store{$tab[0]}{$tab[1]}{"cover"}=$tab[3];
	}
	close $csv;
}


my %store_seeds;
my %store_clustered;
foreach my $genome (sort {$seq_len{$b} <=> $seq_len{$a} || $a cmp $b} keys %seq_len){
	my $c_c='';
	my $id=0;
	my $th_c=$th_cover/100*$seq_len{$genome};
	my $th_match=$th_id*$th_c/100;
	print "$genome - length $seq_len{$genome} - th will be $th_c bp / $th_match matches\n";
	foreach my $seed (sort {$seq_len{$b} <=> $seq_len{$a} || $a cmp $b} keys %store_seeds){
		if (defined($store{$genome}{$seed})){
			if ($store{$genome}{$seed}{"match"}>=$th_match && $store{$genome}{$seed}{"cover"}>=$th_c){
				print "$genome - $seed - ".$store{$genome}{$seed}{"match"}." --- ".$store{$genome}{$seed}{"cover"}."\n";
				my $id_temp=$store{$genome}{$seed}{"match"}/$store{$genome}{$seed}{"cover"}*100;
				if ($id_temp>$id){
					$c_c=$seed;
					$id=$store{$genome}{$seed}{"match"}/$store{$genome}{$seed}{"cover"}*100;
				}
			}
		}
	}
	if ($c_c eq ''){
		$store_seeds{$genome}=1;
		print "$genome is a new seed\n";
	}
	else{
		$store_clustered{$c_c}{$genome}=$id;
		print "$genome is now clustered with $c_c ($store_clustered{$c_c}{$genome})\n";
	}
# 	<STDIN>;
}


my $c_id=0;
open my $s1,">",$out_file;
foreach my $genome (sort {$seq_len{$b} <=> $seq_len{$a}} keys %store_seeds){
	print $s1 ">Cluster_".$c_id."\t$genome\t$seq_len{$genome}\n";
	foreach my $clustered (sort {$seq_len{$b} <=> $seq_len{$a}} keys %{$store_clustered{$genome}}){
		print $s1 "$clustered\t$store_clustered{$genome}{$clustered}\t$seq_len{$clustered}\n";
	}
	$c_id++;
}
close $s1;


if ($tag_fna==0){
	open my $s2,">",$out_file_fasta;
	open my $fa,"<",$fasta_file;
	my $c_c='';
	my $tag=0;
	while(<$fa>){
		chomp($_);
		if ($_=~/^>(\S+)/){
			$c_c=$1;
			$tag=0;
			if ($store_seeds{$c_c}==1){
				print $s2 "$_\n";
				$tag=1;
			}
		}
		elsif($tag==1){
			print $s2 "$_\n";
		}
	}
	close $fa;
	close $s2;
}

sub sort {
	my $file_in=$_[0];
	my $file_out=$_[1];
	my $cmd="gawk -F ' ' '{print \$13\" \"\$0}' $file_in | sort > $file_out";
	print $cmd."\n";
	my $out=`$cmd`;
	print $out."\n";
}
