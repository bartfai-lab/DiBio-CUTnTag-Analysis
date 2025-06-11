#!/bin/bash

#Script to average (2) multiple bedgraph files to merge replicates
#version 2, supporting merging of 3 replicates and checks if bedgraphs are zipped or not before unionbedg

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
      file3=$OPTARG
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

##check if files are compressed and unzip if necessary
if [[ "$file1" == *.gz ]]; then

  echo "unzipping $file1"
  gzip -dk $file1 
  
  file1="${file1%.gz}"
  
  echo "$file1" 
  
else
  echo "not zipped bedgraph for file1"
  
fi

if [[ "$file2" == *.gz ]]; then

  echo "unzipping $file2"
  gzip -dk $file2 
  
  file2="${file2%.gz}"
  
  echo "$file2" 
  
else
  echo "not zipped bedgraph for file1"
  
fi

#### generate common bedgraph with bedtools
### if 3 files are given
if [[ -n "$file3" ]]; then

if [[ "$file3" == *.gz ]]; then

  echo "unzipping $file3"
  gzip -dk $file3 
  
  file3="${file3%.gz}"
  
  echo "$file3" 
  
else
  echo "not zipped bedgraph for file1"
  
fi


  bedtools unionbedg -i ${file1} ${file2} ${file3} -g ${genomeSizes} > ${projPath}${name}.union.bdg


#calculate average and generate new file
  awk '{OFS="\t"; print $1, $2, $3, ($4 + $5 + $6) / 3}' ${projPath}${name}.union.bdg > ${projPath}${name}.avg.bdg

else
  bedtools unionbedg -i ${file1} ${file2} -g ${genomeSizes} > ${projPath}${name}.union.bdg


#calculate average and generate new file
  awk '{OFS="\t"; print $1, $2, $3, ($4 + $5) / 2}' ${projPath}${name}.union.bdg > ${projPath}${name}.avg.bdg
fi

##append header and zip for genome browser
#make begraph header
bdg_header="track type=bedGraph name=${name}.avg visibility=full color=0,0,0 altcolor=105,105,105 autoScale=off viewLimits=0:30 windowingFunction=mean smoothingWindow=6"

echo "$bdg_header" > ${projPath}${name}.avg_withHeader.bdg
cat ${projPath}${name}.avg.bdg >> ${projPath}${name}.avg_withHeader.bdg

gzip -k ${projPath}${name}.avg_withHeader.bdg
