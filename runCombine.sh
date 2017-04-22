#!/bin/bash
#author: Yuansheng Liu

LONG="long.fasta"
OUTPUT="bicolor_long.fasta"
TEMPFOLDER="temp_file"

# make combine
# the first is long read file
# the second is the temporary folder
# the last is output file
./combine ${LONG} ${TEMPFOLDER} ${OUTPUT}

# rm -rf ${TEMPFOLDER}
