# Rabies_gene_annotator
Pipline escrito en bash para la anotación y obtención de los 5 genes del virus de la rabia

# Importante:

Tenga en cuenta que sus secuencias deben estar nombradas algo así; de lo contrario el script fallará al intentar reconocer y formular los ID's: RABV-2025-Rab90.fa

Si lo desea, puede hacer los cambios pertinentes en las lineas de cada ciclo "for", para señalar cuál es el ID de las secuencias: 

# ####################################################################

ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')

inf=$(basename ${assembly} | cut -d '-' -f '2') 

newID=${ID}-${inf}

# #####################################################################

# Preparar la base de datos

En la carpeta "RABV_prot_db" se encuentran las bases de datos curadas de los 5 genes de la rabia. Deberá agregar la ruta donde se encuentren estos archivos y sus archivos generados por makeblast (señalado abajo) a su ~/.bashrc del siguiente modo:

~ export Bx_RABVp_DB_PATH="/home/marcopterix/db/blast_db/RABV" 

Una vez hecho esto, deberá compilar la base de datos de BLAST con:

~ makeblastdb -in Nucleoprot_RABV.fasta -dbtype prot -out ./Nprot_RABV_db




# Dependencias necesarias:

-Deberá tener instalado BLAST+ para poder ejecutar BLASTx (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables) y samtools para poder ejecutar samtools faidx.
