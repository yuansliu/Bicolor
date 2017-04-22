#!/bin/bash
#author: Yuansheng Liu

LONG="long.fasta"
SHORT="ERR022075_smp.fasta"
TEMPFOLDER="temp_file"
CORES=28
STEP=3
rm -rf ${TEMPFOLDER}
mkdir ${TEMPFOLDER}

# make iterative
# the first is long read file
# the second is short read file
# the third is initial length of k-mer
# the fourth is iterative round
# the last is a temporary folder
time ./iterative ${LONG} ${SHORT} 13 ${STEP} ${TEMPFOLDER}
time ./iterative ${LONG} ${SHORT} 15 ${STEP} ${TEMPFOLDER}
time ./iterative ${LONG} ${SHORT} 17 ${STEP} ${TEMPFOLDER}
time ./iterative ${LONG} ${SHORT} 19 ${STEP} ${TEMPFOLDER}

cd ${TEMPFOLDER}

ls | grep "^[0-9]" > name.txt

cp ../muscle ./
chmod 777 muscle

echo "begin muscle... plsease wait it finish"

cat name.txt | xargs -i --max-procs=${CORES} bash -c "./muscle -in {} -out _{} -maxiters 1 -diags -quiet >/dev/null 2>&1"

# ./improved long.fasta ${TEMPFOLDER} bicolor_long.fasta

