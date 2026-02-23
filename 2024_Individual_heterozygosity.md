# WGS: Heterozygosity 

Bamfiles mapped to Adamallo are here:
```
/home/ben/projects/rrg-ben/ben/2024_cliv_allo_WGS/bamz_mapped_to_Adamallo_with_bubbles
```
I ran this script using a list of male and female bam files (so the heteroz is calculated by sex across all individuals within each sex:
```
#!/bin/sh
#SBATCH --job-name=angsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=64:00:00
#SBATCH --mem=512gb
#SBATCH --output=angsd.%J.out
#SBATCH --error=angsd.%J.err
#SBATCH --account=rrg-ben

module load StdEnv/2023 angsd/0.940

angsd -P 1 -setMinChunkSize 50 -minMapQ 30 -minQ 20 -nQueueSize 20 -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/allo.fasta.contigs.fasta.gz -GL 1 -out ${1}_out
realSFS ${1}_out.saf.idx -P 1 -fold 1 > ${1}_out.sfs
realSFS saf2theta ${1}_out.saf.idx -sfs ${1}_out.sfs -outname ${1}_out
thetaStat do_stat ${1}_out.thetas.idx -win 10000 -step 10000  -outnames ${1}_theta.thetasWindow.gz
```


# RNAseq: Calculating pi in windows for each genome
I'm using angsd for this.

bam files are here:
```
/home/ben/scratch/2024_allo_muel_RNAseq
```

First extract a new bam file that has only chr7 for each sample. To do this first make a tab-delimited file with the coordinates of the chromosome that you want to extract; e.g.:
```
Chr7L 1 139837619
```
and name it something meaningful such as `Chr7L_1_139837619.txt`

```
/home/ben/projects/rrg-ben/ben/2024_allo_muel_RNAseq/ben_scripts/2024_samtools_subset_bam.sh
```
```
#!/bin/sh
#SBATCH --job-name=samtools_subset_bam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00 
#SBATCH --mem=2gb
#SBATCH --output=samtools_subset_bam.%J.out 
#SBATCH --error=samtools_subset_bam.%J.err
#SBATCH --account=rrg-ben

module load StdEnv/2023  gcc/12.3 samtools/1.20
samtools view ${1} --region-file ${2} -b > ${1}_${2::-4}.bam
samtools index ${1}_${2::-4}.bam
```

make a file (bamfilename_path.txt) that has the path to the bam file (do this for each one).
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/2023_XT_genomz/ben_scripts/2023_angsd_genomicwindows_pi.sh
```
```
#!/bin/sh
#SBATCH --job-name=angsd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=8gb
#SBATCH --output=angsd.%J.out
#SBATCH --error=angsd.%J.err
#SBATCH --account=def-ben


module load StdEnv/2023 angsd/0.940

#angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta -G
L 1 -out ${1}_angsd_out

#angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2020_XT_v10_refgenome/XENTR_10.0_genome_scafconcat.fasta -GL 1 -out ${1}_out
angsd -bam ${1}.txt -doSaf 1 -anc /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XENLA_10.1_genome.fa -GL 1 -out ${1}_out
realSFS ${1}_out.saf.idx -P 24 -fold 1 > ${1}_out.sfs
realSFS saf2theta ${1}_out.saf.idx -sfs ${1}_out.sfs -outname ${1}_out
thetaStat do_stat ${1}_out.thetas.idx -win 50000 -step 10000  -outnames ${1}_theta.thetasWindow.gz
```

This generated an output file with the suffix ".pestPG" that has the windowed stats.

I need to divide these by the number of sites to get the pi per site...

