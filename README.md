# Bidirectional PBWT (Positional Burrows-Wheeler Transform)

## Introduction
Bidirectional PBWT is an efficient method to find clusters of matches around each variant site and study the changes of matching blocks using forward and reverse PBWT at each variant site at the same time.  The input data for Bidirectional PBWT are phased genotype data (in VCF format). 
## Dependencies
- C++ (at least GCC 5)  
- GNU Make  
- Bash  

## Installation
To install the program clone the repository to a local folder using:

`git clone https://github.com/ZhiGroup/bi-PBWT.git`

Enter the repository folder and compile the program:

`cd bi-PBWT`  
`make`

## Usage Instructions
Type  
`./biPBWT`  
by itself to show the help page.  

An example:  
`./biPBWT --readVcf "example.vcf" --writeTo "output" --length 500000 --gap 0 --checkpoint 100000`  


