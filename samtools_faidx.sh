#!/bin/bash

echo -e "#############################################################################################" "\n"
echo -e ========== Filtrar las regiones de los genes en nucleótidos de secuencias de rabia ========== "\n"
echo -e "\t" ===== Inicio: $(date) ===== "\n"
echo -e "#############################################################################################" "\n"

bash BLASTx_annotate_RABV.sh

#------------------------------------------------------------------------------------
#Definir ubicación de los directorios de los archivos .fa y del directorio de salida
dirfa="/home/bioinfocenasa/Analisis_Corridas/SPAdes/virus" #Directorio donde están los archivos .fa
dirout="/home/bioinfocenasa/Analisis_Corridas/SPAdes/virus/BLASTx_annotate" #Directorio de salida de archivos y los archivos *_start_end_*.tsv generados por BLAST
#------------------------------------------------------------------------------------

cd ${dirout}

for gen in N P M G L; do
echo -e "===== Buscando ${gen} ====="

for i in *_start_end_*.tsv; do
    IDi=$(basename ${i} | cut -d '_' -f '1')
    prot=$(basename ${i} | cut -d '_' -f '4' | cut -d '.' -f '1')
    inic=$(awk '{print $1}' ${i})
    fin=$(awk '{print $2}' ${i})

for a in ${dirfa}/*fa; do
    contig=$(basename ${a} | cut -d '.' -f '1')
    IDa=$(basename ${a} | cut -d '-' -f '3' | cut -d '.' -f '1')
    inf=$(basename ${a} | cut -d '-' -f '2')
    newID=${IDa}-${inf}

#--------------------------
# Identificación del gen N
#--------------------------

case ${gen} in N)
	if [[ ${gen} == ${prot} ]]; then
echo -e "Comparando gen: ${gen} ${prot}"
	if [[ ${IDi} == ${newID} ]]; then
echo -e "Comparando ID: ${IDi} ${newID}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${contig}_${prot}_tmp.fa

sed 's/:.*$/_N/' ${dirout}/${contig}_${prot}_tmp.fa > ${dirout}/${contig}_${prot}.fa

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
echo -e "Comparando gen: ${gen} ${prot}"
        if [[ ${IDi} == ${newID} ]]; then
echo -e "Comparando ID: ${IDi} ${newID}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${contig}_${prot}_tmp.fa

sed 's/:.*$/_P/' ${dirout}/${contig}_${prot}_tmp.fa > ${dirout}/${contig}_${prot}.fa

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
echo -e "Comparando gen: ${gen} ${prot}"
        if [[ ${IDi} == ${newID} ]]; then
echo -e "Comparando ID: ${IDi} ${newID}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${contig}_${prot}_tmp.fa

sed 's/:.*$/_M/' ${dirout}/${contig}_${prot}_tmp.fa > ${dirout}/${contig}_${prot}.fa

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
echo -e "Comparando gen: ${gen} ${prot}"
        if [[ ${IDi} == ${newID} ]]; then
echo -e "Comparando ID: ${IDi} ${newID}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${contig}_${prot}_tmp.fa

sed 's/:.*$/_G/' ${dirout}/${contig}_${prot}_tmp.fa > ${dirout}/${contig}_${prot}.fa

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
echo -e "Comparando gen: ${gen} ${prot}"
        if [[ ${IDi} == ${newID} ]]; then
echo -e "Comparando ID: ${IDi} ${newID}"

samtools faidx ${a} ${contig}:${inic}-${fin} > ${dirout}/${contig}_${prot}_tmp.fa

sed 's/:.*$/_L/' ${dirout}/${contig}_${prot}_tmp.fa > ${dirout}/${contig}_${prot}.fa

	  fi
	fi
     esac
   done
  done
done

rm ${dirfa}/CENASA-MX20-Rab47.fa.fai
rm ${dirout}/*_tmp.fa
rm ${dirout}/*_start_end_*.tsv

mkdir -p ${dirout}/Nucleotidos
mkdir -p ${dirout}/Proteinas

mv ${dirout}/*fa ${dirout}/Nucleotidos
mv ${dirout}/*fasta ${dirout}/Proteinas
