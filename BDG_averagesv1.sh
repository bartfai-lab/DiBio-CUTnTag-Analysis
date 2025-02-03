#!/bin/bash

#Script to average (2) multiple bedgraph files to merge replicates
projPath=""
file1=""
file2=""
file3=""
genomeSizes="/vol/molbio/bartfai/jgockel/genomes/Pf/chrNameLength.txt"
name=""

while getopts ":p:a:b:c:g:n:" opt; do
  case $opt in
    p)
      projPath=$OPTARG
      ;;
    a)
      file1=$OPTARG
      ;;
    b)
      file2=$OPTARG
      ;;
    c)
      file2=$OPTARG
      ;;
    g)
      genomeSizes=$OPTARG
      ;;
    n)
      name=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done



## generate common bedgraph with bedtools
#3rd file not implemented
bedtools unionbedg -i ${file1} ${file2} -g ${genomeSizes} > ${projPath}${name}.union.bdg


#calculate average and generate new file
awk '{OFS="\t"; print $1, $2, $3, ($4 + $5) / 2}' ${projPath}${name}.union.bdg > ${projPath}${name}.avg.bdg


##append header and zip for genome browser
#make begraph header
bdg_header="track type=bedGraph name=${name}.avg visibility=full color=0,0,0 altcolor=105,105,105 autoScale=off viewLimits=0:30 windowingFunction=mean smoothingWindow=6"

echo "$bdg_header" > ${projPath}${name}.avg_withHeader.bdg
cat ${projPath}${name}.avg.bdg >> ${projPath}${name}.avg_withHeader.bdg

gzip -k ${projPath}${name}.avg_withHeader.bdg
