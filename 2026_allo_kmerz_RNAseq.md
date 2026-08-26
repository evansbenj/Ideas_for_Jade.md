# Count kmerz in RNAseq

directory:
```
/home/ben/projects/rrg-ben/ben/2024_allo_muel_RNAseq/fq/allo_fq
```

First get seqs from alignment files that have 29mers spaning a transcript specific SNP. Then extract kmers using this script
```
#!/usr/bin/perl
use strict;
use warnings;

# Check for correct arguments
die "Usage: perl extract_kmers.pl <k_value> <input.fasta>\n" unless @ARGV == 2;

my ($k, $input_file) = @ARGV;

# Open the FASTA file
open(my $fh, '<', $input_file) or die "Could not open file '$input_file' $!";

# Set the record separator to '>' to process FASTA format record-by-record
local $/ = '>';

# Discard the first empty record created by the record separator
my $empty_header = <$fh>;

while (my $record = <$fh>) {
    chomp $record;
    
    # Split the header and the sequence
    my ($header, @seq_lines) = split(/\n/, $record);
    my $sequence = join('', @seq_lines);
    
    # Clean the sequence: remove spaces and newlines
    $sequence =~ s/\s+//g;
    $sequence = uc($sequence);
    
    my $seq_len = length($sequence);
    
    # Skip sequences that are shorter than k
    if ($seq_len < $k) {
        warn "Warning: Sequence '$header' is shorter than k=$k. Skipping.\n";
        next;
    }
    
    # Extract and print the k-mers
    print "--- K-mers for $header ---\n";
    for my $i (0 .. $seq_len - $k) {
        my $kmer = substr($sequence, $i, $k);
        print "$kmer\n";
    }
}

```

# Count kmers in fq files for each sample (it can do F and R reads together:
```
#!/bin/sh
#SBATCH --job-name=makemeryldb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=48:00:00
#SBATCH --mem=128gb
#SBATCH --output=makemeryldb.%J.out
#SBATCH --error=makemeryldb.%J.err
#SBATCH --account=rrg-ben

/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl count ${1}_trim.R[1,2].fq.gz threads=4 memory=128 k=29 output ${1}_trim.R12_meryldb.out
```
# Count reads in target
The target is a fasta file with each kmer from the target seq
```
#!/bin/sh
#SBATCH --job-name=makemeryldb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=makemeryldb.%J.out
#SBATCH --error=makemeryldb.%J.err
#SBATCH --account=rrg-ben


/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl count k=29 allo_XXXW_29merz.fa output allo_XXXW_target_kmers.meryl
```

# Intersect with counts from each sample
```
#!/bin/sh
#SBATCH --job-name=meryl_intersect
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_intersect.%J.out
#SBATCH --error=meryl_intersect.%J.err
#SBATCH --account=rrg-ben


/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl intersect ${1} allo_XXXW_target_kmers.meryl output ${1}_XXX_counts.meryl
```
# print out the counts:
```
#!/bin/sh
#SBATCH --job-name=meryl_print
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=meryl_print.%J.out
#SBATCH --error=meryl_print.%J.err
#SBATCH --account=rrg-ben


/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl print ${1}_XXXW_counts.meryl > ${1}_XXXW_kmers_with_counts.txt
```

This worked perfectly. Only females have kmers from XBW. Yeay!

# Pulling out reads with specific kmers

Not needed for XB but useful for allo:
```
#!/bin/sh
#SBATCH --job-name=makemeryldb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=2:00:00
#SBATCH --mem=128gb
#SBATCH --output=makemeryldb.%J.out
#SBATCH --error=makemeryldb.%J.err
#SBATCH --account=rrg-ben


/home/ben/projects/rrg-ben/ben/2025_bin/meryl/build/bin/meryl-lookup -include \
  -sequence ${1}__trim_R1.fq.gz ${1}__trim_R2.fq.gz \
  -mers allo_XXXW_target_kmers.meryl \
  -output ${1}_XXX_R1.fq.gz ${1}_XXX_R2.fq.gz
```
# dh@@

Expression detected in two tadz only (both male):
```
[ben@l4.nibi allo_fq]$ ls -l *dh*W_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad13_S1_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad14_S2_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben 9086 May 17 16:31 X_allofraseri_tad15_S3_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad16_S4_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad17_S5_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad19_S6_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad20_S7_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben 9164 May 17 16:31 X_allofraseri_tad21_S8_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad23_S9_dh*_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 May 17 16:31 X_allofraseri_tad24_S10_dh*_29mers_with_counts.txt
```
#a@@18
Expression also detected in the two male tadz:
```
[ben@l4.nibi allo_fq]$ ls -l *_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad13_S1_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad14_S2_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben 1120 Aug 25 20:49 X_allofraseri_tad15_S3_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad16_S4_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad17_S5_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad19_S6_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad20_S7_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben 1952 Aug 25 20:49 X_allofraseri_tad21_S8_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad23_S9_aldh18aW_29mers_with_counts.txt
-rw-r----- 1 ben rrg-ben    0 Aug 25 20:49 X_allofraseri_tad24_S10_aldh18aW_29mers_with_counts.txt
```
