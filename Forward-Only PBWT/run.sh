#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Program: PBWT (Positional Burrows-Wheeler Transform)"
	echo ""
	echo "Contact: Shaojie Zhang [shzhang@cs.ucf.edu] or Degui Zhi [degui.zhi@uth.tmc.edu]"
	echo ""
	echo "Usage: ./run.sh [options] parameters"
	echo ""
	echo "Required Parameters:"
	echo -e "\t--readVcf <file>\tVCF file"
	echo -e "\t--map <file>\t\tGenetic Mapping file"
	echo ""
	echo "Optional Parameters:"
	echo -e "\t--writeTo <filename>\tOutput filename and location (parameter can be full file path or just filename) [Default = VCF filename]"
	echo -e "\t--length <integer>\tBlock length (in units of centimorgans (cM)) [Default = 1]"
	echo -e "\t--width <integer>\tNumber of haplotypes in block [Default = 100]"
	exit 1
fi

OPTIONS=r:m:o:l:w:
LONGOPTS=readVcf:,map:,writeTo:,length:,width:

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
eval set -- "$PARSED" 

readVcf="" map="" writeTo="" length=1 width=100 
while true; do
	case "$1" in
		-r|--readVcf)
			readVcf="$2"
			shift 2
			;;
		-m|--map)
			map="$2"
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

if [ "$map" == "" ]; then
	echo "The Genetic Mapping file must be specified with the required option --map <file>"
	exit 1
fi

basename=$(basename $readVcf)
filename="${basename%.*}"
if [ "$writeTo" = "" ]; then
	writeTo="$filename"
fi

echo "Running PBWT..."
./PBWT "$readVcf" "$writeTo" "$length" "$width" "$map"
echo "PBWT Finished."

