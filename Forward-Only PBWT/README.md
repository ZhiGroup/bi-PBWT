## Forward-Only PBWT (Positional Burrows-Wheeler Transform)

## Introduction
The Forward-Only PBWT is an efficient method to find blocks of matches. The Positional Burrows-Wheeler Transform (PBWT) was developed by Richard Durbin as a representation of haplotype data for storing the data and finding matches efficiently among a set of haplotypes. The input data for Forward-Only PBWT is phased genotype data (in VCF format) and a genetic map.

## Dependencies
- C++ (at least GCC 5)  
- GNU Make  
- GNU getopt
- Bash  

## Installation
To install the program clone the repository to a local folder using:

`git clone https://github.com/ZhiGroup/bi-PBWT.git`

Enter the repository folder and compile the program:

`cd "bi-PBWT/Forward-Only PBWT"`  
`make`

## Usage Instructions
Type  
`./run.sh`  
by itself to show the help page.  

|         Option         |                         Parameter                        |                                                    Description                                                    |
|:----------------------:|:--------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------:|
| `-r` or `--readVcf`    | Full or relative file path to input VCF file             | Sets the input VCF file on which to run PBWT.                                                                     |
| `-m` or `--map`    	 | Full or relative file path to Genetic Mapping file       | Sets the Genetic Mapping file.                                                                                    |
| `-o` or `--writeTo`    | Full or relative file path and filename for output files | Sets the location and filename for PBWT output files. The default option is the VCF filename.                     |
| `-l` or `--length`     | Block length (in units of centimorgan (cM))              | Sets the minimum length requirement for blocks. The default value is 1 centimorgan.                               |
| `-w` or `--width`      | Block width                                              | Sets the minimum number of haplotypes required for a block to be reported. The default value is 100 haplotypes.   |

An example:  
`./run.sh --readVcf "example.vcf" --map "example.rmap" --writeTo "output" --length 0.1 --width 500` 

## Genetic Mapping Format
The format of the genetic mapping file must be 2 space-separated fields per line: "site number" "genetic mapping".

## Results
When finished executing, PBWT will generate 1 file with the extensions ".blocks".

The file ".blocks" represents each block on its own line with five space-separated fields "starting genomic location of block" "ending genomic location of block" "starting genetic location of block" "ending genetic location of block" "width of block" followed by space seperated ID's of all the haplotypes in the block. IDs are suffixed with either "-0" or "-1" indicating the first and second haplotype of the individual ID, respectively.
