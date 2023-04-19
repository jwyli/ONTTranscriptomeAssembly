# Scripts to run analysis 

<!-- 
$homedir="/hdd1/home/p21_wyli" 
chmod 755 assemble.sh #if not done
chmod 755 compare.sh #if not done
-->


## Files and Directories:
+ homedir: assemble.sh and compare.sh files here
+ 00.dataset: ${homedir}/00.dataset
+ 01.refgenome: ${homedir}/01.refgenome
+ 02.assembly: ${homedir}/02.assembly
+ 03.sqanti: ${homedir}/03.sqanti
+ 04.gffcompare: ${homedir}/04.gffcompare

## 1. Transcriptome assembly
Uses flair, stringtie (both free and guided by reference annotation) and rnabloom to perform transcriptome assembly from reads.fastq   
The output files are stored in 02.assembly 

### Create environment if not done
``` bash
conda env create -f assemble_env.yml
```

### Activate environment
``` bash
conda activate assemble 
```

### Run script
Change $homedir in assemble.sh if necessary
```bash
cd ${homedir}
./assemble.sh baseForName path/to/reads.fastq path/to/refGenome.fasta path/to/refAnnotation.gtf
```
### Output
Output from each sample will be stored in a separate dir baseForName
+ baseForName.flair.collapse.isoforms.fa
+ baseForName.flair.collapse.isoforms.gtf
+ baseForName.stringtie.free.assembly.fa
+ baseForName.stringtie.free.assembly.gtf
+ baseForName.stringtie.guided.assembly.fa
+ baseForName.stringtie.guided.assembly.gtf
+ baseForName.rb.transcripts.fa (no gtf for de novo assembly)

### Run RATTLE separately
RATTLE takes very long and preferable to run it in a separate wrapper script.

### Activate environment
```bash
conda activate rattle
export LD_LIBRARY_PATH=~/6010/source/RATTLE/spoa/build/lib
```
### Run script
```bash
cd ${homedir}
./runrattle.sh baseForName path/to/reads.fasta
```

## 2. Assessment of assembly quality
Uses sqanti and gffcompare to compare the assembled transcriptome (assembly.gtf | transcripts.fa) with reference transcriptome (refAnnotation.gtf)

### Create environment if not done
Easiest way is create an env to install SQANTI first then install gffcompare in the same env  
(Follow https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation) 


### Activate environment
``` bash
conda activate cmp
```

### Run script
change $homedir and $sqantiPath in compare.sh if necessary
```bash
cd ${homedir}/
./compare.sh baseForName path/to/refGenome.fasta path/to/refAnnotation.gtf
```

### Output
The output files are stored in 03.sqanti and 04.gffcompare under dir baseForName  
+ sqanti uses .fa files for analysis and outputs results in frFa dir  
(sqanti also creates a .renamed.fasta file in 02.assembly for its own use, where the transcript description line is changed)
+ gffcompare only uses available .gtf files for analysis (ie. only stringtie and flair assemblies)

## 3. Summarizing sqanti results
run sqanti_summarize.sh after going to the directory with the sqanti output
```bash
${homedir}/sqanti_summarize.sh
```
Output: tsv with number of transcripts in different sqanti categories for each assembly

run missing.sh after going to directory with the sqanti output 
```bash
${homedir}/missing.sh /path/to/refAnnotation.gtf
```
Output:
+ refGenes.tmp: list of genes in the reference (ground truth) set
+ refTranscripts.tmp: list of transcripts in the reference (ground truth) set
+ $assembler.genes.tmp: list of genes associated with the transcripts assembled
+ $assembler.transcripts.tmp: list of transcripts assembled
+ $assembler.genes.diff: list of missing genes in the assembly (compared with reference set)
+ $assembler.transcripts.diff: list of missing transcripts in the assembly (compared with reference set)

###### Temp link for results generated 
+ Dataset SRR6058583 Oxford Nanopore direct RNA of SIRV mix E2: https://www.ncbi.nlm.nih.gov/sra/?term=SRR6058583  
+ Reference genome and annotations: https://www.lexogen.com/wp-content/uploads/2021/06/SIRV_Set1_Norm_Sequences_20210507.zip
+ Selected files from generated results: https://www.dropbox.com/scl/fo/tu1au2ogd7o9fd289xd08/h?dl=0&rlkey=heir5i03ne0ayfuqd00ykn52g



