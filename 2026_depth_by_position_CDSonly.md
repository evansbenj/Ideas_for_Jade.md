# Depth by position CDS only

Because repetitive seqs mess things up with depth, I'm interested in trying to get depth for individual positions but only within CDS.

To accomplish this I used this XL file as a query against the boum genome:
```
/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_boum_genome_assembly/XL_CDS_only_nospaces.fasta
```

Then I parsed this output file to make a bed file with this script:
```perl
#!/usr/bin/perl
use warnings;
use strict;

# This program will read in a blast output file and output a bed file that
# has the coordinates of all the matches that are at least 90% identical
# and also at least 90% of the length of the query 

# I am doing this because I will use this bed file as input to samtools
# to get depth of each pygm WGS sample mapped to the boum genome

# My hope is to then identify putative coding regions where the coverage in all females
# is high and the coverage in all males is low

my $inputfile = $ARGV[0];
my $outputfile1 = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";



my @temp;
my $percentage = 0.9;
my @query1;
my @query2;
my @query3;
my %hash=();
my $key;

while ( my $line = <DATAINPUT>) {
	@temp = split("\t",$line);
	print "hi $temp[0]\n";
	@query1 = split("_",$temp[0]);
	@query2 = split(":",$query1[$#query1]);
	@query3 = split("-",$query2[$#query2]);
	# check if the hit is at least $percentage long and $percentage identical
	if($query3[1]>$query3[0]){
		if(($temp[2]>=($percentage*100))&&($temp[3]>= ($query3[1]-$query3[0])*$percentage)){
			# this is a good hit so save it
			$hash{$temp[0]."_".$temp[1]."_".$temp[8]}[0]= $temp[1]; # boum tig
			$hash{$temp[0]."_".$temp[1]."_".$temp[8]}[1]= $temp[8]; # boum start
			$hash{$temp[0]."_".$temp[1]."_".$temp[8]}[2]= $temp[9]; # boum end
		}	
	}
	else{
		if(($temp[2]>=($percentage*100))&&($temp[3]>= ($query3[0]-$query3[1])*$percentage)){
			# this is a good hit so save it
			$hash{$temp[0]."_".$temp[1]."_".$temp[8]}[0]= $temp[1]; # boum tig
			$hash{$temp[0]."_".$temp[1]."_".$temp[8]}[1]= $temp[8]; # boum start
			$hash{$temp[0]."_".$temp[1]."_".$temp[8]}[2]= $temp[9]; # boum end
		}
	}	
}	
close DATAINPUT;


foreach $key (sort(keys %hash)) {
	if($hash{$key}[1] < $hash{$key}[2]){
		print OUTFILE $hash{$key}[0],"\t",$hash{$key}[1],"\t",$hash{$key}[2],"\n";
	}else{
		print OUTFILE $hash{$key}[0],"\t",$hash{$key}[2],"\t",$hash{$key}[1],"\n";
	}	
}
close OUTFILE;
```
