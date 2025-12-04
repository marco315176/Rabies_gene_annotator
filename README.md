# Rabies_gene_annotator
Pipline escrito en bash para la anotación y obtención de los 5 genes individuales del virus de la rabia en formato fasta  


# Importante:

Tenga en cuenta que sus secuencias deben estar nombradas algo así; de lo contrario el script fallará al intentar reconocer y formular los ID's: **RABV-2025-Rab90.fa**

Por ahora, los contigs de su genoma completo de rabia deben llamarse igual que su archivo fasta; de modo que, siguiendo el ejemplo anterior, su contig debe estár identificado como: **>RABV-2025-Rab90**

Si lo desea, puede hacer los cambios pertinentes en las lineas de cada ciclo "for" en los dos scripts aquí proporcionados, para señalar cuál es el ID de las secuencias: 

```
#Cambiando:
ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')

inf=$(basename ${assembly} | cut -d '-' -f '2') 

newID=${ID}-${inf}

#Por algo como:
ID=$(basename ${assembly} | cut -d '-' -f '1') #Si su archivo fasta se nombra algo como: Rab90-assembly.fa
```

# Dependencias necesarias:

***Deberá tener instalado:   
-> BLAST+ para poder ejecutar BLASTx (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables)  
-> samtools para poder ejecutar samtools faidx (https://www.htslib.org/) y  
-> seqkit (https://bioinf.shenwei.me/seqkit/download/)***

# Preparar la base de datos

Para descargar las bases de datos deberás tener instalado previamente BLAST+ y seqkit. Una vez que los tenga instalados, se debe ejecutar el script RGA_db_dwl.sh del siguiente modo:

```
bash RGA_db_dwl.sh
```

Esto descargará las bases de datos en la carpeta $HOME/db/RGA. Posteriormente, deberá agregar la ruta donde se encuentren estos archivos generados a su ~/.bashrc del siguiente modo:

```
export Bx_RABV_RGA_PATH="$HOME/db/RGA"
```


# Ejecutar el pipline

Una vez que tenga todos los requisitos y haya configurado las rutas de entrada y salida en los scripts .sh, ejecute el siguiente script, mismo que ejecutará tanto **BLASTx_annotate_RABV.sh** como la identificación con **samtools faidx**: 

```
bash rabies_gene_annotator.sh
```

# Archivos de salida

Al final de la ejecución del pipline, en la carpeta que usted marque como ${dirout}, encontrará dos carpetas: Nucleotidos y Proteinas, en las cuales encontrará las secuencias de nucleotidos y aminoácidos, respectivamente; así como también encontrará el archivo **Annotation_all.tsv** donde estará la información de la anotación de sus secuencias.



