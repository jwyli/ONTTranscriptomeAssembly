#!/bin/bash

homedir="/hdd1/home/p21_wyli/6010/"
baseForName=$1 
readsFq=$2 #full path

exec > >(tee "${baseForName}.rattle.out") 2>&1 # save all subsequent output
set -x #echo on

starttime=`date +%s`

cd ${homedir}/02.assembly
mkdir ${baseForName}
cd ${baseForName}

/hdd1/home/p21_wyli/6010/source/RATTLE/rattle cluster -i ${readsFq} -t 48 --iso --rna 
/hdd1/home/p21_wyli/6010/source/RATTLE/rattle correct -i ${readsFq} -c "${homedir}/02.assembly/${baseForName}/clusters.out" -t 48 
/hdd1/home/p21_wyli/6010/source/RATTLE/rattle polish -i consensi.fq -t 48 --rna

endtime=`date +%s`
runtime=$((endtime-starttime))
echo "finished RATTLE analysis for ${baseForName} in ${runtime} seconds."

mv clusters.out "${baseForName}.rattle.clusters.out"
mv corrected.fq "${baseForName}.rattle.corrected.fq"
mv uncorrected.fq "${baseForName}.rattle.uncorrected.fq"
mv consensi.fq "${baseForName}.rattle.consesnsi.fq"
mv transcriptome.fq "${baseForName}.rattle.transcriptome.fq"
