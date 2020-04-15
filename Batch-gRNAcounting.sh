#! /bin/bash


conda activate mypy3

for file in *R2*.fastq.gz
do
    recom_file=${file/R2/R_2_recom}
    ../DesktopMacBookPro/lab/0331_lib_distribution/seqkit seq -r -p $file | gzip -c  > $recom_file

    gunzip $recom_file

    awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 150) {print header, seq, qheader, qseq}}' < ${recom_file/fastq.gz/fastq} > ${recom_file/fastq.gz/filtered.fastq}


   # sleep 3
    python ../DesktopMacBookPro/lab/0331_lib_distribution/count_spacers.py -f ${recom_file/fastq.gz/filtered.fastq} -no-g
   # sleep 3
    mv library_count.csv ${recom_file/fastq.gz/csv}
   # sleep 3
done