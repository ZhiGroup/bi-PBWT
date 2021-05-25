# Bidirectional PBWT (Positional Burrows-Wheeler Transform)

## Introduction
Bidirectional PBWT is an efficient method to find clusters of matches around each variant site and study the changes of matching blocks using forward and reverse PBWT at each variant site at the same time. The Positional Burrows-Wheeler Transform (PBWT) was developed by Richard Durbin as a representation of haplotype data for storing the data and finding matches efficiently among a set of haplotypes. The input data for Bidirectional PBWT is phased genotype data (in VCF format).

## Dependencies
- C++ (at least GCC 5)  
- GNU Make  
- GNU getopt
- Bash  
- [Optional] Boost C++ Libraries (required for memory mapped files)

## Installation
To install the program clone the repository to a local folder using:

`git clone https://github.com/ZhiGroup/bi-PBWT.git`

Enter the repository folder and compile the program:

`cd bi-PBWT`  
`make`

[Optional] Compile with `make MemoryMap` if using memory mapped files.

## Usage Instructions
Type  
`./biPBWT.sh`  
by itself to show the help page.  

|         Option         |                         Parameter                        |                                                    Description                                                    |
|:----------------------:|:--------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------:|
| `-r` or `--readVcf`    | Full or relative file path to input VCF file             | Sets the input VCF file on which to run biPBWT.                                                                   |
| `-o` or `--writeTo`            | Full or relative file path and filename for output files | Sets the location and filename for biPBWT output files. The default option is the VCF filename.                   |
| `-l` or `--length`     | Block length (in units of base pairs)                    | Sets the minimum length requirement for blocks. The default value is 500,000 base pairs.                          |
| `-w` or `--width`      | Block width                                              | Sets the minimum number of haplotypes required for a block to be reported. The default value is 100 haplotypes.   |
| `-g` or `--gap`        | Gap size (in units of sites)                             | Sets the size of the gap in the block. The default value is 1 site.                                               |
| `-c` or `--checkpoint` | Checkpoint interval `n`                                  | biPBWT will print a checkpoint message to console every `n` sites. The default value is 100,000.                  |
| `-s` or `--sites`              | `None`                                                   | Change the units of distance for length to sites. The default is base pairs.                                      |

An example:  
`./biPBWT.sh --readVcf "example.vcf" --writeTo "output" --length 500000 --width 100 --gap 0 --checkpoint 100000`  

## Results
When finished executing, biPBWT will generate 3 files with the extensions ".clusters", ".IDs", and ".MI".

The file ".clusters" represents each cluster on its own line with seven space seperated fields "site k" "genomic location of site k" "forward length (in sites)" "reverse length (in sites)" "starting genomic location of block" "ending genomic location of block" "width of block".

The file ".IDs" will have the same number of lines as the ".clusters" file with each line representing the ID's of all the haplotypes in the cluster on the same line in the ".cluster" file. IDs are suffixed with either "-1" or "-2" indicating the first and second haplotype of the individual ID, respectively.

The file ".MI" has two space seperated values on each line representing "genomic location" "mutual information value".

Note that biPBWT generates 3 intermediate files with the extensions ".rpbwt", ".sites", and ".meta" during execution that are deleted upon completion. The ".rpbwt" file stores the reverse positional prefix array and reverse divergence array in binary format (each value takes 4 bytes) and takes up around twice the disk space of the VCF file. The ".sites" file stores the genomic position of each site on its own line. The ".meta" files stores two space seperated values M and N representing the number of haplotypes and the number of sites in the VCF file.

## Citation
If you found our work useful in your research, please consider citing the following paper:

