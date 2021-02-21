#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Program: Bidirectional PBWT (Positional Burrows-Wheeler Transform)"
	echo ""
	echo "Contact: Shaojie Zhang [shzhang@cs.ucf.edu] or Degui Zhi [degui.zhi@uth.tmc.edu]"
	echo ""
	echo "Usage: ./biPBWT.sh [options] parameters"
	echo ""
	echo "Required Parameters:"
	echo -e "\t--readVcf <file>\tVCF file"
	echo ""
	echo "Optional Parameters:"
	echo -e "\t--writeTo <filename>\tOutput filename and location (parameter can be full file path or just filename) [Default = VCF filename]"
	echo -e "\t--length <integer>\tBlock length (in units of base pairs) [Default = 500000]"
	echo -e "\t--width <integer>\tNumber of haplotypes in block [Default = 100]"
	echo -e "\t--gap <integer>\t\tGap Size (site) [Default = 1]"
	echo -e "\t--checkpoint <integer>\tConsole output every n sites [Default = 100000]"
	echo -e "\t--sites\t\t\tChange units of distances for length to sites [Default = base pairs]"
	exit 1
fi

OPTIONS=c:r:o:l:w:g:s
LONGOPTS=checkpoint:,readVcf:,writeTo:,length:,width:,gap:,sites

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
eval set -- "$PARSED" 

checkpoint=100000 readVcf="" writeTo="" length=500000 width=100 gap=1 sites=0
while true; do
	case "$1" in
		-c|--checkpoint)
			checkpoint="$2"
			shift 2
			;;
		-r|--readVcf)
			readVcf="$2"
			shift 2
			;;
		-o|--writeTo)
			writeTo="$2"
			shift 2
			;;
		-l|--length)
			length="$2"
			shift 2
			;;
		-w|--width)
			width="$2"
			shift 2
			;;
		-g|--gap)
			gap="$2"
			shift 2
			;;
		-s|--sites)
			sites=1
			shift
			;;
		--)
			shift
			break
			;;
	esac
done

if [ "$readVcf" == "" ]; then
	echo "The VCF input file must be specified with the required option --readVcf <file>"
	exit 1
fi

basename=$(basename $readVcf)
filename=$(echo $basename | cut -f1 -d ".")
if [ "$writeTo" = "" ]; then
	writeTo="$filename"
fi

echo "Running rPBWT..."
./rPBWT "$readVcf" "$writeTo" "$checkpoint"
echo "Running PBWT..."
./PBWT "$readVcf" "$writeTo" "$checkpoint" "$length" "$width" "$gap" "$sites"
rm "${writeTo}.rpbwt" "${writeTo}.sites" "${writeTo}.meta"
echo "biPBWT Finished."

