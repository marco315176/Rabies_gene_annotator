# Rabies_gene_annotator

Pipline escrito en bash para la anotación y obtención de los 5 genes individuales del virus de la rabia en nucleótidos y aminoácidos.


# Importante:

Este pipeline funciona bien con secuencias completas del virus de la rabia (RABV), por lo que es importante mencionar que si usted quisiera anotar secuencias fragmentadas, samtools podría llegar a fallar al intentar extraer las secuencias de sus genes. Esto se debe a que, al realizar la anotación de los genes con secuencias de proteinas como referencia, para incluir el codón de paro, se busca tres bases después de lo indicado por BLASTx; de modo que si el contig en el que se anotó el gen es más pequeño, no obtendrá la secuencia de este.

```
#working...
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

Esto descargará las bases de datos en la carpeta $HOME/db/RGA. Si gusta puede agregar la ruta donde se encuentren estos archivos generados a su ~/.bashrc del siguiente modo:

```
export Bx_RABV_RGA_PATH="$HOME/db/RGA"
```


# Ejecutar el pipline

Una vez que tenga todos los requisitos, deberá cambiar los permisos de sus scripts para que puedan ejecutarse, de este modo:

```
chmod +x RGA_*
chmod +x rabies_gene_annotator.sh
```
Y para ejecutar el pipeline: 

```
bash rabies_gene_annotator.sh
```

# Archivos de salida

Al final de la ejecución del pipline, en la carpeta que usted marque como ${dirout}, encontrará dos carpetas: Nucleotidos y Proteinas, en las cuales encontrará las secuencias de nucleotidos y aminoácidos, respectivamente; así como también encontrará el archivo **Annotation_all.tsv** donde estará la información de la anotación de sus secuencias.



