# Mapping RNAseq data to allofraseri genome assembly

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

