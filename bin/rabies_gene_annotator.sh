#!/usr/bin/env bash

echo -e "\n""\n""\033[43m========== RABV gene annotator ==========\033[m"
echo -e "\033[43m====== $(date) =====\033[m""\n"

#---------
# Options
#---------

usage () {
echo ""
echo -e "\033[4;33m========== Bash script designed to annotate and extract amino acid and nucleotide sequences of the five RABV genes. ==========\033[0m"
echo ""
echo "Pipeline developed in the massive sequencing and bioinformatics area of CENASA, SENASICA"
echo ""
echo -e "\033[4;33mThis script executes first and second steps of this Rabies Gene Annotator pipeline.\033[0m"
echo ""
echo "Options:"
echo "Usage: $0 -f FASTA file PATH -o OUTDIR PATH"
echo " -h print help "
echo " -f FASTA file directory "
echo " -o OUTPUT directory "
echo " -p PATH to BLAST database. If you downloades the database by running the RGA_db_dwl.sh script, the path to your database is: $HOME/db/RGA "
echo "";
	}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

while getopts ":hf:o:p:" opt; do
      case ${opt} in
h)
  usage; exit
;;
f)
  dirfa=${OPTARG}
;;
o)
  dirout=${OPTARG}
;;
p)
  dirdb=${OPTARG}
;;
:)
  echo -e "\033[0;33mOption -${OPTARG} requires an argument.\033[0m"
exit
;;
\?)
  echo -e "\033[0;31mInvalid option -${OPTARG}.\033[0m"
exit 1
;;
	esac
     done


if [[ ${dirfa} == ${dirout} ]]; then
    echo -e "\033[0;31mError: -f and -o cannot be the same PATH.\033[0m"
   exit 1
fi

if [[ -z ${dirfa} ]]; then
    echo -e "\033[0;33mError: Option -f is necesary.\033[0m"
    exit 1
fi

if [[ -z ${dirout} ]]; then
    echo -e "\033[0;33mError: Option -o is necesary.\033[0m"
    exit 1
fi

if [[ -z ${dirdb} ]]; then
    echo -e "\033[0;33mError: Option -p is necesary.\033[0m"
    exit 1
fi

#------------------------------
# Run BLASTx to annotate genes
#------------------------------

RGA_BLASTx_annotate.sh -f ${dirfa} \
                       -o ${dirout} \
                       -p ${dirdb}

#
# Run samtools to extract the nucleotide sequence of the genes
#

RGA_samtools_faidx.sh -f ${dirfa} \
                      -o ${dirout}

echo -e "\033[3;32m###############################################################\033[0m"
echo -e "\033[3;32m========== RABV genome annotation with RGA completed ==========\033[0m"
echo -e "\033[3;32m================= $(date) ================\033[0m"
echo -e "\033[3;32m###############################################################\033[0m"
