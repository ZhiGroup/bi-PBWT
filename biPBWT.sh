#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Program: Bidirectional PBWT (Positional Burrows-Wheeler Transform)"
	echo ""
	echo "Contact: Shaojie Zhang [shzhang@cs.ucf.edu] or Degui Zhi [degui.zhi@uth.tmc.edu]"
	echo ""
	echo "Usage: ./biPBWT.sh [options]"
	echo ""
	echo "Required Parameters:"
	echo -e "\t--readVcf <file>\tVCF file"
	echo ""
	echo "Optional Parameters:"
	echo -e "\t--writeTo <filename>\tOutput filename and location"
	echo -e "\t--length <integer>\tBlock Length (bp) [Default = 500000]"
	echo -e "\t--gap <integer>\t\tGap Size (site) [Default = 0]"
	echo -e "\t--checkpoint <integer>\tConsole output every n sites [Default = 100000]"
	exit 1
fi

OPTIONS=c:r:w:l:g:
LONGOPTS=checkpoint:,readVcf:,writeTo:,length:,gap:

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
eval set -- "$PARSED" 

checkpoint=100000 readVcf="" writeTo="" length=50 gap=0
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
		-w|--writeTo)
			writeTo="$2"
			shift 2
			;;
		-l|--length)
			length="$2"
			shift 2
			;;
		-g|--gap)
			gap="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
	esac
done

basename=$(basename $readVcf)
filename=$(echo $basename | cut -f1 -d ".")
if [ "$writeTo" = "" ]; then
	writeTo="$filename"
fi

echo "Running rPBWT..."
./rPBWT "$readVcf" "$writeTo" "$checkpoint"
echo "Running PBWT..."
./PBWT "$readVcf" "$writeTo" "$checkpoint" "$length" "$gap"
rm "${writeTo}.rpbwt" "${writeTo}.sites" "${writeTo}.meta"
echo "biPBWT Finished."
