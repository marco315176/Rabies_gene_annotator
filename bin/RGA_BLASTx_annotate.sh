#!/usr/bin/env bash

echo -e "\n""\n""\033[42m=============== Annotating ==============\033[m"
echo -e "\033[42m====== $(date) =====\033[m""\n"

#---------
# Options
#---------

usage () {
echo ""
echo -e "\033[4;33m===== Bash script designed to annotate and extract amino acid and nucleotide sequences of the five RABV genes. =====\033[0m"
echo ""
echo "Pipeline developed in the massive sequencing and bioinformatics area of CENASA, SENASICA"
echo ""
echo -e "\033[4;33mThis is the first step: Annotation with BLASTx and obtaining amino acid sequences.\033[0m"
echo ""
echo "Options:"
echo "Usage: $0 -f FASTA file PATH -o OUTDIR PATH -p BLAST DB PATH"
echo " -h print help "
echo " -f FASTA file directory "
echo " -o OUTPUT directory "
echo " -p PATH to BLAST database. If you downloaded the database by running the RGA_db_dwl.sh script, the path to your database is: $HOME/db/RGA "
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

#
 if [[ ${dirfa} == ${dirout} ]]; then
     echo -e "\033[0;31mError: -f and -o cannot be the same PATH.\033[0m"
     exit 1
fi

#

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

mkdir -p ${dirout}
cd ${dirfa}

#--------------------------------------
#Definir DB de proteinas de los 5 genes
#--------------------------------------

dbN="Nprot_RABV_db"
dbP="Pprot_RABV_db"
dbM="Mtxprot_RABV_db"
dbG="Gprot_RABV_db"
dbL="Lprot_RABV_db"

#-------------------------------------
#Loop para encontrar el gen N de RABV
#-------------------------------------

echo -e "\n\033[1;36m========== Searching for gene N ==========\033[0m\n"

for assembly in *.fa *.fasta *.fna; do
    ID=$(basename ${assembly} | cut -d '.' -f '1')

echo -e "##### ${ID} #####"

blastx -query ${assembly} \
       -db ${dirdb}/${dbN} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 \
       -out ${dirout}/${ID}_N_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_N""\n"$11}' ${dirout}/${ID}_N_blastx.tsv \
> ${dirout}/${ID}_N_prot.fna

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""N"}' ${dirout}/${ID}_N_blastx.tsv \
> ${dirout}/${ID}_N_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen P de RABV
#--------------------------------------

echo -e "\n\033[1;36m========== Searching for gene P ==========\033[0m\n"

for assembly in *.fa *.fasta *.fna; do
    ID=$(basename ${assembly} | cut -d '.' -f '1')

echo -e "##### ${ID} #####"

blastx -query ${assembly} \
       -db ${dirdb}/${dbP} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 \
       -out ${dirout}/${ID}_P_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_P""\n"$11}' ${dirout}/${ID}_P_blastx.tsv \
> ${dirout}/${ID}_P_prot.fna

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""P"}' ${dirout}/${ID}_P_blastx.tsv \
> ${dirout}/${ID}_P_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen M de RABV
#--------------------------------------

echo -e "\n\033[1;36m========== Searching for gene M ==========\033[0m\n"

for assembly in *.fa *.fasta *.fna; do
    ID=$(basename ${assembly} | cut -d '.' -f '1')

echo -e "##### ${ID} #####"

blastx -query ${assembly} \
       -db ${dirdb}/${dbM} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 \
       -out ${dirout}/${ID}_M_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_M""\n"$11}' ${dirout}/${ID}_M_blastx.tsv \
> ${dirout}/${ID}_M_prot.fna

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""M"}' ${dirout}/${ID}_M_blastx.tsv \
> ${dirout}/${ID}_M_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen G de RABV
#--------------------------------------

echo -e "\n\033[1;36m========== Searching for gene G ==========\033[0m\n"

for assembly in *.fa *.fasta *.fna; do
    ID=$(basename ${assembly} | cut -d '.' -f '1')

echo -e "##### ${ID} #####"

blastx -query ${assembly} \
       -db ${dirdb}/${dbG} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 \
       -out ${dirout}/${ID}_G_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_G""\n"$11}' ${dirout}/${ID}_G_blastx.tsv \
> ${dirout}/${ID}_G_prot.fna

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""G"}' ${dirout}/${ID}_G_blastx.tsv \
> ${dirout}/${ID}_G_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen L de RABV
#--------------------------------------

echo -e "\n\033[1;36m========== Searching for gene L ==========\033[0m\n"

for assembly in *.fa *.fasta *.fna; do
    ID=$(basename ${assembly} | cut -d '.' -f '1')

echo -e "##### ${ID} #####"

blastx -query ${assembly} \
       -db ${dirdb}/${dbL} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 \
       -out ${dirout}/${ID}_L_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_L""\n"$11}' ${dirout}/${ID}_L_blastx.tsv \
> ${dirout}/${ID}_L_prot.fna

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""L"}' ${dirout}/${ID}_L_blastx.tsv \
> ${dirout}/${ID}_L_info.tsv

done


#--------------------------------------------------------------
# Concatenar archivos resultantes y eliminar los que no se usan
#--------------------------------------------------------------

rm ${dirout}/*_blastx.tsv

cd ${dirout}

for n in *N_info.tsv; do
    IDn=$(basename ${n} | cut -d '_' -f '1')
    p=${n/_N_/_P_}
    IDp=$(basename ${p} | cut -d '_' -f '1')
    m=${n/_N_/_M_}
    IDm=$(basename ${m} | cut -d '_' -f '1')
    g=${n/_N_/_G_}
    IDg=$(basename ${g} | cut -d '_' -f '1')
    l=${n/_N_/_L_}
    IDl=$(basename ${l} | cut -d '_' -f '1')

	if [[ ${IDn} != ${IDp} && ${IDm} && ${IDg} && ${IDl} ]]; then
continue
	else

#Concatenar los archivos de información del alineamiento en uno solo por ID
cat ${n} ${p} ${m} ${g} ${l} > ${dirout}/${IDn}_info_align.tsv

sed -i '1i ID\tReferencia\tIdentidad\tAling_long\tInicio\tFin\te-value\tbitscore\tGen' ${dirout}/${IDn}_info_align.tsv

	fi

done

#----------------------------------------
# Concatenar archivos .fna de proteinas
#----------------------------------------

for n in *_N_prot.fna; do
    IDn=$(basename ${n} | cut -d '_' -f '1')
    p=${n/_N_/_P_}
    IDp=$(basename ${p} | cut -d '_' -f '1')
    m=${n/_N_/_M_}
    IDm=$(basename ${m} | cut -d '_' -f '1')
    g=${n/_N_/_G_}
    IDg=$(basename ${g} | cut -d '_' -f '1')
    l=${n/_N_/_L_}
    IDl=$(basename ${l} | cut -d '_' -f '1')

 if [[ ${IDn} != ${IDp} && ${IDm} && ${IDg} && ${IDl} ]]; then
continue
        else

#Concatenar los archivos de archivos fasta en uno solo por ID
cat ${n} ${p} ${m} ${g} ${l} > ${dirout}/${IDn}_prot_RABV.fna

	fi
done

rm *_prot.fna


echo -e "\033[5;32m######################################################################\033[0m"
echo -e "\033[5;32m========== Annotation of gene regions with BLASTx completed ==========\033[0m"
echo -e "\033[5;32m==================== $(date) ====================\033[0m"
echo -e "\033[5;32m######################################################################\033[0m"


