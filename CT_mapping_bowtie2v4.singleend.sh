#!/usr/bin/env bash
#this is version 4 for single end ChIP data
#define variables
projPath=""
filename=""
genome=""
duplicates=TRUE
mapQuality=""
bigWig=TRUE
exChr="PFC10_API_IRAB M76611"
MIT="PFC10_API_IRAB"
API="M76611"


# parse command line options and arguments
while getopts ":n:g:p:d:b:x:q:" opt; do
  case $opt in
    n)
      filename=$OPTARG
      ;;
    g)
      genome=$OPTARG
      ;;
    p)
      projPath=$OPTARG
      ;;
    d)
      duplicates=$OPTARG
      ;;
    b)
      bigWig=$OPTARG
      ;;
    x)
      exChr=$OPTARG
      ;;
    q)
      mapQuality=$OPTARG
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

#defines filename with -n, genome location with -g, -q filter quality

# shift the parsed options so that the remaining arguments are at the front of the list
shift $((OPTIND-1))


#generate necessary directories to store alignment files
mkdir -p ${projPath}/${filename}
mkdir -p ${projPath}/${filename}/alignment/sam/
mkdir -p ${projPath}/${filename}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/${filename}/alignment/bam
mkdir -p ${projPath}/${filename}/alignment/bed
mkdir -p ${projPath}/${filename}/alignment/bedgraph
mkdir -p ${projPath}/${filename}/alignment/bigWig



#add indexing step potentially
############################PLACEHOLDER##############################################

#map with bowtie2 for samples build correctly under projPath and filename

#check if bowtie2 has run before
if [ -e ${projPath}/${filename}/alignment/sam/${filename}_bowtie2.sam ]; then
    echo "Data has been mapped, skip bowtie2 alignment"
else
    echo "performing bowtie2 mapping"

#Run bowtie2
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 8 -x ${genome} -U ${projPath}/fastq/${filename}.fastq.gz -S ${projPath}/${filename}/alignment/sam/${filename}_bowtie2.sam &> ${projPath}/${filename}/alignment/sam/bowtie2_summary/${filename}_bowtie2.txt

echo "done with bowtie2 mapping"

fi

projPath=${projPath}/${filename}

###########################################################################################

filename2=${filename}

#generate filter report
echo "Start of the Run: $(date)" >> $projPath/alignment/${filename2}.filter_report.txt
echo -e "$(samtools view -c $projPath/alignment/sam/${filename}_bowtie2.sorted.sam)\treads before filtering" >> $projPath/alignment/${filename2}.filter_report.txt

#######################################################################################

#remove duplicates
if [[ "$duplicates" == TRUE ]]; then
  echo "perform duplicate removal with picard"
  
  #generate directory
  mkdir -p $projPath/alignment/removeDuplicate/picard_summary
  
  ## Sort by coordinate
  picard SortSam -I $projPath/alignment/sam/${filename}_bowtie2.sam -O $projPath/alignment/sam/${filename}_bowtie2.sorted.sam -SORT_ORDER coordinate


  ## mark duplicates
  picard MarkDuplicates -I $projPath/alignment/sam/${filename}_bowtie2.sorted.sam -O $projPath/alignment/removeDuplicate/${filename}_bowtie2.sorted.dupMarked.sam -METRICS_FILE $projPath/alignment/removeDuplicate/picard_summary/${filename}_picard.dupMark.txt

  ## remove duplicates
  picard MarkDuplicates I=$projPath/alignment/sam/${filename}_bowtie2.sorted.sam O=${projPath}/alignment/sam/${filename}.rmDup_bowtie2.sorted.sam REMOVE_DUPLICATES=true METRICS_FILE=$projPath/alignment/removeDuplicate/picard_summary/${filename}_picard.rmDup.txt
  
 
  #add line to filter report
	echo -e "$(samtools view -c ${projPath}/alignment/sam/${filename}.rmDup_bowtie2.sorted.sam)\treads left after removing duplicates" >> $projPath/alignment/${filename2}.filter_report.txt
  
  filename=${filename}.rmDup

  
else
  echo "skip duplicate removal"

if [ -e $projPath/alignment/sam/${filename}_bowtie2.sorted.sam ]; then
    echo "sam already sorted"
else
    
echo "sorting sam with picard"
  ## Sort by coordinate
  picard SortSam -I $projPath/alignment/sam/${filename}_bowtie2.sam -O $projPath/alignment/sam/${filename}_bowtie2.sorted.sam -SORT_ORDER coordinate
  
fi
fi

#################################################################################################

#filter reads by mapping quality

if [ -e ${projPath}/alignment/sam/${filename}_bowtie2.sorted.qualityScore$mapQuality.sam ]; then
    echo "bams have already been filtered by map quality"
    
    #add line to filter report
	echo -e "$(samtools view -c ${projPath}/alignment/sam/${filename}_bowtie2.sorted.qualityScore$mapQuality.sam)\treads left after filtering for $mapQuality map quality" >> $projPath/alignment/${filename2}.filter_report.txt
 
 
else
    
echo "Filter by $mapQuality map quality"

samtools view -h -q $mapQuality ${projPath}/alignment/sam/${filename}_bowtie2.sorted.sam >${projPath}/alignment/sam/${filename}_bowtie2.sorted.qualityScore$mapQuality.sam

#add line to filter report
	echo -e "$(samtools view -c ${projPath}/alignment/sam/${filename}_bowtie2.sorted.qualityScore$mapQuality.sam)\treads left after filtering for $mapQuality map quality" >> $projPath/alignment/${filename2}.filter_report.txt
  
fi

#####################################################################################################

#file conversion into bams

if [ -e $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.bam ]; then
    echo "bams have already been generated, skipping ahead"
else
    
echo "generate bams"



## Filter and keep the mapped read pairs on min quality scores
samtools view -bS -F 0x04 "${projPath}/alignment/sam/${filename}_bowtie2.sorted.qualityScore$mapQuality.sam" > $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.bam

##sort bam
samtools sort $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.bam > $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.bam

fi

########################################################################################################

#works until here

########################################################################################################

#remove apicoplast and mitrochondrial reads

#extract header from bam file
samtools view -H $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.bam -> $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.header.tmp

#remove Mitochondria and Apicoplast
#generate sam file that is able to be filtered
samtools view $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.bam > $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp 
								
#filter out chroms 
		echo "filtering out ${MIT}"
		
			# filter out out chrom with grep -w -v chrom
			grep -w -v ${MIT} $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp > $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.sam.tmp
			cp $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.sam.tmp $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp
			
			echo "done removing ${MIT}"
			#add line to filter report
			echo -e "$(wc -l $chrom $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp | cut -f1 -d' ')\treads left after removing ${MIT}" >> $projPath/alignment/${filename2}.filter_report.txt


echo "filtering out ${API}"
			
			# filter out out chrom with grep -w -v chrom
			grep -w -v ${API} $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp > $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.sam.tmp
			cp $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.sam.tmp $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp
			
			echo "done removing ${API}"
			#add line to filter report
			echo -e "$(wc -l $chrom $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp | cut -f1 -d' ')\treads left after removing ${API}" >> $projPath/alignment/${filename2}.filter_report.txt



		#concatinate header and sam file
		cat $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.header.tmp $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp |  samtools view -bS -> $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bam

#remove temporary files
rm $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.header.tmp
rm $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.sam.tmp 
rm $projPath/alignment/bam/${filename}_bowtie2.mapped.qualityScore$mapQuality.sam.tmp


################################################################################
#works until here
################################################################################

#generate bigwigs
if [ "$bigWig" == TRUE ]; then

  if [ -e $projPath/alignment/bigWig/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.RPMK.bw ]; then
    echo "bigWigs have already been generated, skipping ahead"
  else
    
  echo "bigWigs being generated"

  samtools index $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bam
  
  bamCoverage --bam $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bam -o $projPath/alignment/bigWig/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.RPMK.bw --normalizeUsing RPKM
  fi 
  
else

  echo "bigWig production skipped"

fi

#################################################################################

##compute scaling factor
#lbrary size
librarysize=$(samtools view -c  $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bam)
echo "library size is: $librarysize"

#to per M
scale_per_M=$( echo 1000000/${librarysize} | bc -l )



#make begraph header
bdg_header="track type=bedGraph name=${filename}.Q${mapQuality}.perM.bdg description=${filename}.Q${mapQuality}.perM.bdg visibility=full color=0,0,0 altcolor=105,105,105 autoScale=off viewLimits=0:50 windowingFunction=mean smoothingWindow=4"

#check for bedgraphs
if [ -e $projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bedgraph ]; then
    echo "bedgraphs have already been generated, skipping ahead"
else
    
echo "generate bedgraphs"

# Calculate genome coverage
if bedtools genomecov -scale ${scale_per_M} -bg -ibam "$projPath/alignment/bam/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bam" > "$projPath/alignment/bedgraph/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bedgraph"; then
    echo "Genome coverage calculation successful."
else
    echo "Error: Genome coverage calculation failed."
    exit 1
fi

# Prepend BEDGraph header
if { echo -e "${bdg_header}"; cat "$projPath/alignment/bedgraph/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit.bedgraph"; } > "$projPath/alignment/bedgraph/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit_with_header.bedgraph"; then
    echo "BEDGraph header prepended successfully."
else
    echo "Error: Failed to prepend BEDGraph header."
    exit 1
fi

# Compress the resulting file
if gzip -k "$projPath/alignment/bedgraph/${filename}_bowtie2.mapped.sorted.qualityScore$mapQuality.noApiMit_with_header.bedgraph"; then
    echo "File compression successful."
else
    echo "Error: File compression failed."
    exit 1
fi

echo "Done producing bedgraph."

fi

