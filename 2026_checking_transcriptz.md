# Checking transcriptz

Based on gene content, I think these 11 are Y-linked: 
* tig00006146
* tig00005690
* tig00005092
* tig00005113
* tig00005103
* tig00005106
* tig00006042
* tig00006035
* tig00011671
* tig00009702
* tig00009703

And these are X-linked:
* tig00010978
* tig00009681

To identify transcripts with a perfect match to these contigs, I made blastdbs from each one:

```
tig00005092	1	567460
tig00005103	1	548174
tig00005106	1	526359
tig00005113	1	918335
tig00005690	1	2208206
tig00006035	1	130300
tig00006042	1	544463
tig00006146	1	1416245
tig00009702	1	359031
tig00009703	1	386102
tig00011671	1	1575154
```

```
bedtools getfasta -fi allo.fasta.contigs.fasta -bed Y_tigs.bed -fo Y_tigs.fa
```
```
makeblastdb -in Y_tigs.fa -dbtype nucl -out Y_tigs.fa_blastable
```
Before blasting the RNAseq transcripts I added an underscore so the length is part of the name:
```
sed -i 's/ len=/_len=/g' allo_trinity_assembly_all_batches.Trinity.fasta
```
Then, the blast output is prefiltered to include only 100% matches like this:
```
grep '      100.000 ' allo_denovotranscripts_to_Y_tigs > expressed_Y_trnascriptz.out
```

```
blastn -query allo_trinity_assembly_all_batches.Trinity.fasta -db Y_tigs.fa_blastable -outfmt 6 -out allo_denovotranscripts_to_Y_tigs
```

Parse the results using this perl script
```perl
#!/usr/bin/perl
use warnings;
use strict;

# This program reads in a blast output from a RNAseq transcriptome assembly to a WGS genome assembly.
# Before blasting the RNAseq transcripts I added an underscore so the length is part of the name:
# sed -i 's/ len=/_len=/g' allo_trinity_assembly_all_batches.Trinity.fasta

# Then, the blast output is prefiltered to include only 100% matches like this:
# grep '      100.000 ' allo_denovotranscripts_to_Y_tigs > expressed_Y_trnascriptz.out

# For each transcript, I will calculate the proportion of the length with 100% match
# I will also flag transcripts with multiple matches to different contigs

# here is the format of the blast output

# Column headers:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

#  1.	 qseqid	 query (e.g., gene) sequence id
#  2.	 sseqid	 subject (e.g., reference genome) sequence id
#  3.	 pident	 percentage of identical matches
#  4.	 length	 alignment length
#  5.	 mismatch	 number of mismatches
#  6.	 gapopen	 number of gap openings
#  7.	 qstart	 start of alignment in query
#  8.	 qend	 end of alignment in query
#  9.	 sstart	 start of alignment in subject
#  10.	 send	 end of alignment in subject
#  11.	 evalue	 expect value
#  12.	 bitscore	 bit score

# ungapped alignment is col 8 - col 7 +1 -col6 (perl col7-col6 +1 -col5)
# do not use col4 because the alignment length includes gaps in the query and ref
# mismatches of nongapped seqs is column 5 (perl col4)

my $outputfile = "expressed_Y_trnascriptz_summary.txt";
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile  $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

my $blast_results = $ARGV[0]; # blast results outfm6 expressed_Y_trnascriptz.txt

unless (open DATAINPUT, $blast_results) {
	print "Can not find the data input file!\n";
	exit;
}


my @temp;
my @temp2;
my %blast_results;  # this will summarize the results.
					# key will be transcript ID
					# value will be total ungapped alignment [0]
					# ungapped mismatches [1]
					# SL or not [2] (1, 0, respectively)
					# earliest [3] and latest position [4] in transcript with a match to the genome
					# the transcript positions will be updated for each section
my $key;
my $previous_gene="";
my $gene_temp;
my $gene;
my $previous_chr="";
my $allowable_overlap=20;
my $firstline=0;


print OUTFILE "transcript\tlength\tmatch_length\tproportion\tgapornot\n";
print "transcript\tlength\tmatch_length\tproportion\tgapornot\n";

while ( my $line = <DATAINPUT>) {
	@temp = split("\t",$line);
	$gene_temp = $temp[0];
	@temp2 = split("_len=",$gene_temp);
	$gene = $temp2[0];
	# check if we are on the same gene 
	# ranges are $temp[6] and $temp[7] and $previous_start and $previous_stop
	if (($gene eq $previous_gene)||($firstline == 0)){ # we are on the same RNA transcript
		# check if there is overlap and 
		# if the next bit is on the same chr
		# with any of the previous ranges
		if(($temp[1] eq $previous_chr)&&($firstline == 1)){ # we are on the same contig and this is not the first line
		# check if the hit overlaps with the previous match
		# only record it if there is not an overlap
			if($temp[6] > $blast_results{$gene}[3]){ # ok this is a new region
				$blast_results{$gene}[1] += $temp[3]; # add this is the length of the match
				$blast_results{$gene}[2] = $temp[6]; # this is the start of the match
				$blast_results{$gene}[3] = $temp[7]; # this is the end of the match
				if($temp[5] > 0){
					$blast_results{$gene}[4] = 1; # this means that there is a gap
				}
				else{
					$blast_results{$gene}[4] = 0; # this means that there is not a gap
				}	
			}
			$previous_gene=$gene;
			$previous_chr=$temp[1];	
		}
		else{ # this is a new contig or the first line; set the flag
			if($firstline == 1){ # this is not the first line but it is another match to a transcript we saw previously
				$blast_results{$gene}[5] = 1; # this transcript hits or spans multiple contigs	
				if($temp[6] > $blast_results{$gene}[3]){ # ok this is a new region
					$blast_results{$gene}[1] += $temp[3]; # add this is the length of the match
					$blast_results{$gene}[2] = $temp[6]; # this is the start of the match
					$blast_results{$gene}[3] = $temp[7]; # this is the end of the match
					if($temp[5] > 0){
						$blast_results{$gene}[4] = 1; # this means that there is a gap
					}
					else{
						$blast_results{$gene}[4] = 0; # this means that there is not a gap
					}
				}		
			}
			else{ # this is the first line
				$blast_results{$gene}[0] = $temp2[1]; # this is the length of the RNA transcript
				$blast_results{$gene}[1] += $temp[3]; # add this is the length of the match
				$blast_results{$gene}[2] = $temp[6]; # this is the start of the match
				$blast_results{$gene}[3] = $temp[7]; # this is the end of the match
				if($temp[5] > 0){
					$blast_results{$gene}[4] = 1; # this means that there is a gap
				}
				else{
					$blast_results{$gene}[4] = 0; # this means that there is not a gap
				}
				$blast_results{$gene}[5] = 0; # this transcript does not (yet)hit or spans multiple contigs	
				print $temp2[0]," $firstline hello\n";			
			}
			$previous_gene=$gene;
			$previous_chr=$temp[1];	
		}
		$firstline=1;
	}	
	
	else{
		if($firstline == 0){
			$firstline = 1;
		}
		else{
			# print out the info from the previous gene
			print OUTFILE $gene," ",$blast_results{$previous_gene}[0],"\t",$blast_results{$previous_gene}[0],"\t", sprintf("%.2f",$blast_results{$previous_gene}[1]/$blast_results{$previous_gene}[0]),"\t",$blast_results{$previous_gene}[4],"\n";
		}
		$blast_results{$gene}[0] = $temp2[1]; # this is the length of the RNA transcript
		$blast_results{$gene}[1] = $temp[3]; # this is the length of the match
		$blast_results{$gene}[2] = $temp[6]; # this is the start of the match
		$blast_results{$gene}[3] = $temp[7]; # this is the end of the match
		if($temp[5] > 0){
			$blast_results{$gene}[4] = 1; # this means that there is a gap
		}
		else{
			$blast_results{$gene}[4] = 0; # this means that there is not a gap
		}
		$blast_results{$gene}[5] = 0;
		$previous_gene=$gene;
		$previous_chr=$temp[1];	
	}
}	

close DATAINPUT;

```

Do the same for X chr and compare
