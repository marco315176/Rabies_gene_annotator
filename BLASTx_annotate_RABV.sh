#!/bin/bash

echo -e "################################################################" "\n"
echo -e ========== Anotación de regiones de genes con BLASTx ========== "\n"
echo -e "\t"              ===== Inicio: $(date) ===== "\n"
echo -e "################################################################" "\n"

#----------------------------------------------------------------------------------------------------------------
#Diseñado para extraer secuencias de proteinas de una base de datos de proteinas comparado con uno de nucleótidos
#----------------------------------------------------------------------------------------------------------------

cd /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus

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

echo -e ========== Buscando gen N ========== "\n"

for assembly in *.fa; do
    ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')
    inf=$(basename ${assembly} | cut -d '-' -f '2')
    newID=${ID}-${inf}

echo -e "##### ${newID} #####"

blastx -query ${assembly} \
       -db $Bx_RABVp_DB_PATH/${dbN} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -evalue 1e-10 \
       -out /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_N_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_N""\n"$11}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_N_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_N_prot.fasta

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""N"}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_N_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_N_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen P de RABV
#--------------------------------------

echo -e ========== Buscando gen P ========== "\n"

for assembly in *.fa; do
    ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')
    inf=$(basename ${assembly} | cut -d '-' -f '2')
    newID=${ID}-${inf}

echo -e "##### ${newID} #####"

blastx -query ${assembly} \
       -db $Bx_RABVp_DB_PATH/${dbP} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -evalue 1e-10 \
       -out /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_P_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_P""\n"$11}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_P_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_P_prot.fasta

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""P"}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_P_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_P_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen M de RABV
#--------------------------------------

echo -e ========== Buscando gen M ========== "\n"

for assembly in *.fa; do
    ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')
    inf=$(basename ${assembly} | cut -d '-' -f '2')
    newID=${ID}-${inf}

echo -e "##### ${newID} #####"

blastx -query ${assembly} \
       -db $Bx_RABVp_DB_PATH/${dbM} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -evalue 1e-10 \
       -out /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_M_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_M""\n"$11}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_M_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_M_prot.fasta

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""M"}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_M_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_M_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen G de RABV
#--------------------------------------

echo -e ========== Buscando gen G ========== "\n"

for assembly in *.fa; do
    ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')
    inf=$(basename ${assembly} | cut -d '-' -f '2')
    newID=${ID}-${inf}

echo -e "##### ${newID} #####"

blastx -query ${assembly} \
       -db $Bx_RABVp_DB_PATH/${dbG} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -evalue 1e-10 \
       -out /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_G_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_G""\n"$11}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_G_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_G_prot.fasta

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""G"}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_G_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_G_info.tsv

done

#--------------------------------------
# Loop para encontrar el gen L de RABV
#--------------------------------------

echo -e ========== Buscando gen L ========== "\n"

for assembly in *.fa; do
    ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')
    inf=$(basename ${assembly} | cut -d '-' -f '2')
    newID=${ID}-${inf}

echo -e "##### ${newID} #####"

blastx -query ${assembly} \
       -db $Bx_RABVp_DB_PATH/${dbL} \
       -max_target_seqs 1 -max_hsps 1 -culling_limit 1 -evalue 1e-10 \
       -out /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_L_blastx.tsv \
       -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend evalue bitscore qseq'

#Extraer proteina en formato FASTA
awk '{print ">"$1"_L""\n"$11}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_L_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_L_prot.fasta

#Extraer información del alineamiento
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7"\t"($8+3)"\t"$9"\t"$10"\t""L"}' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_L_blastx.tsv \
> /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${newID}_L_info.tsv

done


#--------------------------------------------------------------
# Concatenar archivos resultantes y eliminar los que no se usan
#--------------------------------------------------------------

rm /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/*_blastx.tsv

cd /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate

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

#Generar los archivos de inicio y fin del gen que serán usados por samtools faidx
awk '{print $5"\t"$6}' ${n} > ${IDn}_start_end_N.tsv
awk '{print $5"\t"$6}' ${p} > ${IDp}_start_end_P.tsv
awk '{print $5"\t"$6}' ${m} > ${IDm}_start_end_M.tsv
awk '{print $5"\t"$6}' ${g} > ${IDg}_start_end_G.tsv
awk '{print $5"\t"$6}' ${l} > ${IDl}_start_end_L.tsv

	if [[ ${IDn} != ${IDp} && ${IDm} && ${IDg} && ${IDl} ]]; then
continue
	else

#Concatenar los archivos de información del alineamiento en uno solo por ID
cat ${n} ${p} ${m} ${g} ${l} > /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${IDn}_info_align.tsv

sed -i '1i ID\tReferencia\tIdentidad\tAling_long\tInicio\tFin\te-value\tbitscore\tGen' /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${IDn}_info_align.tsv

	fi

done

rm /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/*_info.tsv

#----------------------------------------
# Concatenar archivos .fasta de proteinas
#----------------------------------------

for n in *_N_prot.fasta; do
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
cat ${n} ${p} ${m} ${g} ${l} > /home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate/${IDn}_prot_RABV.fasta

	fi
done

rm *_prot.fasta
