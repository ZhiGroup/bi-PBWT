#!/bin/bash

if [ $# -eq 0 ]; then
	echo "Program: FastRecomb"
fi

vcf="" writeTo="" length=1 width=500 rMap=""
while getopts ":i:o:L:W:m:d:" opt; do
	case ${opt} in
		i)
			vcf=$OPTARG
			;;
		o)	
			writeTo=$OPTARG
			;;
		L)
			length=$OPTARG
			;;
		W)
			width=$OPTARG
			;;
		m)
			rMap=$OPTARG
			;;
		\?)
			echo "Invalid option: $OPTARG"
			;;
		:)
			echo "Invalid option: $OPTARG requries an argument"
			;;
	esac
done
shift $((OPTIND -1))

basename=$(basename $vcf)
filename=$(echo $basename | cut -f1 -d ".")
if [ "$writeTo" = "" ]; then
	writeTo="$filename"
fi

echo "Running PBWT..."
./PBWT "$vcf" "$writeTo" "$length" "$width" "$rMap"
echo "FastRecomb Finished."

