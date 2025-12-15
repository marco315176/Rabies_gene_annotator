#!/usr/bin/env bash

echo -e "\n""\n""\033[42m========== Filter gene regions into rabies sequence nucleotides after annotating ==========\033[m"
echo -e "\033[42m=============================== $(date) ==============================\033[m""\n"

#---------
# Options
#---------

usage () {
echo ""
echo -e "\033[4;33m===== Bash script designed to annotate and extract amino acid and nucleotide sequences of the five RABV genes. =====\033[0m"
echo -e ""
echo "Pipeline developed in the massive sequencing and bioinformatics area of CENASA, SENASICA"
echo ""
echo -e "\033[4;33mThis is the second step: Obtaining nucleotide sequences with samtools faidx, based on the annotation with BLASTx.\033[0m"
echo ""
echo "Options:"
echo "Usage: $0 -f FASTA file PATH -o OUTDIR PATH"
echo " -h print help"
echo " -f FASTA file directory. The same as indicated in the BLASTx_annotate_RABV.sh script."
echo " -o OUTPUT directory. The same as indicated in the BLASTx_annotate_RABV.sh script, again."
echo "";
	}

if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

while getopts ":hf:o:" opt; do
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


###

cd ${dirout}

for gen in N P M G L; do
echo -e "\n\033[1;36m========== Searching ${gen} ==========\033[0m\n"


for i in *info.tsv; do
    IDi=$(basename ${i} | cut -d '_' -f '1')
    prot=$(basename ${i} | cut -d '_' -f '2')
    contig=$(awk '{print $1}' ${i})
    inic=$(awk '{print $5}' ${i})
    fin=$(awk '{print $6}' ${i})

for a in ${dirfa}/*.f*; do
    [[ $a == *.f*.* ]] && continue
    ID=$(basename ${a} | cut -d '.' -f '1')

#--------------------------
# Identificación del gen N
#--------------------------

case ${gen} in N)
        if [[ ${gen} == ${prot} ]]; then

        if [[ ${IDi} == ${ID} ]]; then

echo -e "Comparando ID: ${IDi} ${ID}"
echo -e "Comparando gen: ${gen} ${prot}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${ID}_${prot}_tmp.fa

sed 's/:.*$/_N/' ${dirout}/${ID}_${prot}_tmp.fa > ${dirout}/${ID}_${prot}.fa

else
	continue
   fi
fi
;;

#--------------------------
# Identificación del gen P
#--------------------------

               P)
        if [[ ${gen} == ${prot} ]]; then

        if [[ ${IDi} == ${ID} ]]; then

echo -e "Comparando ID: ${IDi} ${ID}"
echo -e "Comparando gen: ${gen} ${prot}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${ID}_${prot}_tmp.fa

sed 's/:.*$/_P/' ${dirout}/${ID}_${prot}_tmp.fa > ${dirout}/${ID}_${prot}.fa

else
        continue
   fi
fi
;;

#--------------------------
# Identificación del gen M
#--------------------------

               M)
        if [[ ${gen} == ${prot} ]]; then

        if [[ ${IDi} == ${ID} ]]; then

echo -e "Comparando ID: ${IDi} ${ID}"
echo -e "Comparando gen: ${gen} ${prot}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${ID}_${prot}_tmp.fa

sed 's/:.*$/_M/' ${dirout}/${ID}_${prot}_tmp.fa > ${dirout}/${ID}_${prot}.fa

else
        continue
   fi
fi
;;

#--------------------------
# Identificación del gen G
#--------------------------

               G)
        if [[ ${gen} == ${prot} ]]; then

        if [[ ${IDi} == ${ID} ]]; then

echo -e "Comparando ID: ${IDi} ${ID}"
echo -e "Comparando gen: ${gen} ${prot}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${ID}_${prot}_tmp.fa

sed 's/:.*$/_G/' ${dirout}/${ID}_${prot}_tmp.fa > ${dirout}/${ID}_${prot}.fa

else
        continue
   fi
fi
;;

#--------------------------
# Identificación del gen L
#--------------------------

               L)
        if [[ ${gen} == ${prot} ]]; then

        if [[ ${IDi} == ${ID} ]]; then

echo -e "Comparando ID: ${IDi} ${ID}"
echo -e "Comparando gen: ${gen} ${prot}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${ID}_${prot}_tmp.fa

sed 's/:.*$/_L/' ${dirout}/${ID}_${prot}_tmp.fa > ${dirout}/${ID}_${prot}.fa

	  fi
	fi
     esac
   done
  done
 done

rm ${dirfa}/*.fai
rm ${dirout}/*_tmp.fa
rm ${dirout}/*_info.tsv


mkdir -p ${dirout}/Nucleotidos
mkdir -p ${dirout}/Proteinas

mv ${dirout}/*fa ${dirout}/Nucleotidos
mv ${dirout}/*fna ${dirout}/Proteinas

for f in ${dirout}/*info_align.tsv; do
    ename=$(basename ${f} | cut -d '_' -f '1')

echo -e "\n########## ${ename} ########## \n$(cat ${f})"
	done >> ${dirout}/Annotation_all.tsv

rm ${dirout}/*_info_align.tsv


echo -e "\033[3;32m###############################################################\033[0m"
echo -e "\033[3;32m========== Obtaining nucleotide sequences completed. ==========\033[0m"
echo -e "\033[3;32m================= $(date) ================\033[0m"
echo -e "\033[3;32m###############################################################\033[0m"
