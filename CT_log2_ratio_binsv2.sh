#!/usr/bin/env bash

#define variables
background=""
cutntag=""
outputdir=""
binsize=""
outname=""

# parse command line options and arguments
while getopts ":b:c:s:o:n:" opt; do
  case $opt in
    b)
      background=$OPTARG
      ;;
    c)
      cutntag=$OPTARG
      ;;
    s)
      binsize=$OPTARG
      ;;
    o)
      outputdir=$OPTARG
      ;;
    n)
      outname=$OPTARG
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


# shift the parsed options so that the remaining arguments are at the front of the list
shift $((OPTIND-1))

###################################################

#make directories
mkdir ${outputdir}/multiBigWigsummary
mkdir ${outputdir}/multiBigWigsummary/${binsize}binsize

#run multibigwigsummary

if [ -e ${outputdir}/multiBigWigsummary/${binsize}binsize/${binsize}${outname}bin.enrichmentscore.tab ]; then
    echo "multiBigWigsummary has been run, skip ahead"
else
    echo "performing multiBigWigsummary"

multiBigwigSummary bins -b ${cutntag}  ${background} --binSize=${binsize} -o ${outputdir}/multiBigWigsummary/${binsize}binsize/${binsize}bin.results.npz --outRawCounts ${outputdir}/multiBigWigsummary/${binsize}binsize/${binsize}${outname}bin.raw.enrichmentscore.tab


echo "MuliBigWig summary completed"

fi

##################################################

##adapt enrichmentscore file

#remove mito and apicoplast regions even tho they are 0 
awk 'BEGIN {OFS = FS = "\t"} $1 != "M76611" && $1 != "PFC10_API_IRAB" {print}' ${outputdir}/multiBigWigsummary/${binsize}binsize/${binsize}${outname}bin.raw.enrichmentscore.tab > ${outputdir}/multiBigWigsummary/${binsize}binsize/${outname}${binsize}bin.clean.enrichmentscore.tab 


#add 1 pseudocount to all to avoid dividing by 0 and calculate log2 ratio
awk 'BEGIN {OFS = FS = "\t"} {print $1, $2, $3, $4+1, $5+1, log(($4+1)/($5+1))/log(2)}' ${outputdir}/multiBigWigsummary/${binsize}binsize/${outname}${binsize}bin.clean.enrichmentscore.tab > ${outputdir}/multiBigWigsummary/${binsize}binsize/${binsize}${outname}bin.enrichmentscore.pseudocount.tab

##generate bedgraph
echo "generate bedgraph"


awk -F'\t' '{print $1 "\t" $2 "\t" $3 "\t" $6}' ${outputdir}/multiBigWigsummary/${binsize}binsize/${binsize}${outname}bin.enrichmentscore.pseudocount.tab > ${outputdir}/multiBigWigsummary/${binsize}binsize/${outname}.bdg


echo -e "track type=bedGraph name=${outname}.${binsize}binsize description=${outname}.${binsize}binsize visibility=full color=0,0,0 altcolor=105,105,105 autoScale=off viewLimits=-2:4 windowingFunction=mean smoothingWindow=4" > tmpfile
cat ${outputdir}/multiBigWigsummary/${binsize}binsize/${outname}.bdg >> tmpfile
mv tmpfile ${outputdir}/multiBigWigsummary/${binsize}binsize/${outname}.bdg

gzip -k ${outputdir}/multiBigWigsummary/${binsize}binsize/${outname}.bdg
