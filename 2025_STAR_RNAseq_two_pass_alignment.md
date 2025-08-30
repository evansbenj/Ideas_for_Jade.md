# Mapping RNAseq data to allofraseri genome assembly

Directory:
```
/home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/ben_scripts
```

We lack an annotation file for the allofraseri genome but we can still align the RNAseq data to this genome using STAR: https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/

# Or as described by AI (gasp):
```
Step-by-step Guide:
1. Generate the genome index (1st Pass):
Use the --runMode genomeGenerate command with your FASTA file. 
Do NOT include a --sjdbGTFfile at this stage, as you do not have one. 
2. Run the 1st Pass Alignment for all samples:
Run STAR for each of your samples using the index created in the previous step. 
Specify the --runMode alignReads parameter. 
Collect the SJ.out.tab file from each sample. 
3. Combine and filter junction files:
Concatenate all the individual SJ.out.tab files into a single file. 
Apply any necessary filtering (e.g., remove junctions that are too short, or too few reads), though for multiple samples you may not need to filter much. 
4. Generate the 2nd Pass Genome Index:
Create a new index using the genomeGenerate mode, but this time provide the combined SJ.out.tab file using the --sjdbFileChrStartEnd parameter. 
5. Run the 2nd Pass Alignment:
Perform the 2nd pass alignment for all your samples using this new, updated genome index. 
This pass will use the comprehensive set of novel junctions identified in the first pass, improving mapping accuracy for these junctions.
```

# First pass index:
```
#!/bin/sh
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=256gb
#SBATCH --output=star_index.%J.out
#SBATCH --error=star_index.%J.err
#SBATCH --account=rrg-ben

module load star/2.7.11b

STAR --runMode genomeGenerate --genomeDir /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/without_bubbles/ --genomeFastaFiles /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/without_bubbles/allo.fasta.contigs_nobubbles.fasta --runThreadN 8 --limitGenomeGenerateRAM=124544990592
```
# First pass alignment
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=64gb
#SBATCH --output=STAR_map.%J.out
#SBATCH --error=STAR_map.%J.err
#SBATCH --account=def-ben


module load star/2.7.11b

STAR --runMode alignReads \
     --genomeDir /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/with_bubbles/ \
     --runThreadN 6 \
     --readFilesIn ${1} ${2} \
     --outFileNamePrefix ${3} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard \
     --readFilesCommand zcat
```
# Concatenate SJ tab filez
```
cat *tab > allo_all_SJ.out.tab
```
# Second pass index:
```
#!/bin/sh
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=256gb
#SBATCH --output=star_index.%J.out
#SBATCH --error=star_index.%J.err
#SBATCH --account=rrg-ben

module load star/2.7.11b

STAR --runMode genomeGenerate --genomeDir /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/with_bubbles/ --genomeFastaFiles /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/with_bubbles/allo.fasta.contigs.fasta --runThreadN 8 --sjdbFileChrStartEnd allo_all_SJ.out.tab --limitGenomeGenerateRAM=124544990592
```
# Second pass alignment
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=64gb
#SBATCH --output=STAR_map.%J.out
#SBATCH --error=STAR_map.%J.err
#SBATCH --account=def-ben


module load star/2.7.11b

STAR --runMode alignReads \
     --genomeDir /home/ben/projects/rrg-ben/ben/2025_allo_PacBio_assembly/Adam_allo_genome_assembly/with_bubbles/ \
     --runThreadN 6 \
     --readFilesIn ${1} ${2} \
     --outFileNamePrefix ${3} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard \
     --readFilesCommand zcat
```
