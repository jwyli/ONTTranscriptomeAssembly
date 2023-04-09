#!/bin/bash

#on conda activate msqanti should have set $PYTHONPATH for cDNA_cupcake

homedir="/hdd1/home/p21_wyli/6010"
sqantiPath="/hdd1/home/p21_wyli/6010/source/SQANTI3"

baseForName=$1 #baseForName=SRR6058584
refGenome=$2 #refGenome=/hdd1/home/p21_wyli/6010/01.refgenome/SIRV_set1.fasta
refAnnotation=$3 #refAnnotation=/hdd1/home/p21_wyli/6010/01.refgenome/SIRV_set1.gtf
#dataPath=$(dirname "${readsFq}") ; dataFile=$(basename "${readsFq}")

flairGtf="${homedir}/02.assembly/${baseForName}/${baseForName}.flair.collapse.isoforms.gtf" #line 36
stFreeGtf="${homedir}/02.assembly/${baseForName}/${baseForName}.stringtie.free.assembly.gtf" #line 73
stGuidedGtf="${homedir}/02.assembly/${baseForName}/${baseForName}.stringtie.guided.assembly.gtf" #line 74

flairFa="${homedir}/02.assembly/${baseForName}/${baseForName}.flair.collapse.isoforms.fa"
stFreeFa="${homedir}/02.assembly/${baseForName}/${baseForName}.stringtie.free.assembly.fa"
stGuidedFa="${homedir}/02.assembly/${baseForName}/${baseForName}.stringtie.guided.assembly.fa" 
rbFa="${homedir}/02.assembly/${baseForName}/${baseForName}.rb.transcripts.fa"
rattleFa="${homedir}/02.assembly/${baseForName}/${baseForName}.rattle.transcriptome.fq"

exec > >(tee "${baseForName}.cmp.out") 2>&1 # save all subsequent output
set -x #echo on

#perform gffcompare with assembled transcriptomes (flair and stringtie only)
cd ${homedir}/04.gffcompare

mkdir $baseForName
cd $baseForName 

echo "Performing gffcompare for flair and stringtie (free and guided)"

#usage: gffcompare [options]* {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}
gffcompare -r "$refAnnotation" -o "${baseForName}.flair.gffcmp" "$flairGtf" 
gffcompare -r "$refAnnotation" -o "${baseForName}.stringtie.free.gffcmp" "$stFreeGtf"
gffcompare -r "$refAnnotation" -o "${baseForName}.stringtie.guided.gffcmp" "$stGuidedGtf"

#consider -R in some samples, but not sirv

#tidy up
cd "${homedir}/02.assembly/${baseForName}" 
mv *gffcmp* "${homedir}/04.gffcompare/${baseForName}" #move the gffcmp files from 02.assembly back to 04. gffcompare

#perform sqanti3 for assembled transcriptomes 
cd "${homedir}/03.sqanti"
mkdir ${baseForName}
cd ${baseForName}
mkdir frFa

#use fasta files
python "${sqantiPath}/sqanti3_qc.py" $flairFa $refAnnotation $refGenome \
-d "${homedir}/03.sqanti/${baseForName}/frFa" --skipORF --report both --fasta --force_id_ignore

python "${sqantiPath}/sqanti3_qc.py" $stFreeFa $refAnnotation $refGenome \
-d "${homedir}/03.sqanti/${baseForName}/frFa" --skipORF --report both --fasta --force_id_ignore

python "${sqantiPath}/sqanti3_qc.py" $stGuidedFa $refAnnotation $refGenome \
-d "${homedir}/03.sqanti/${baseForName}/frFa" --skipORF --report both --fasta --force_id_ignore

python "${sqantiPath}/sqanti3_qc.py" $rbFa $refAnnotation $refGenome \
-d "${homedir}/03.sqanti/${baseForName}/frFa" --skipORF --report both --fasta --force_id_ignore

python "${sqantiPath}/sqanti3_qc.py" $rattleFa $refAnnotation $refGenome \
-d "${homedir}/03.sqanti/${baseForName}/frFa" --skipORF --report both --fasta --force_id_ignore