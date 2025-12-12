# Rabies_gene_annotator

Pipline escrito en bash para la anotación y obtención de los 5 genes individuales del virus de la rabia en nucleótidos y aminoácidos.


# Importante:

Este pipeline funciona bien con secuencias completas del virus de la rabia (RABV), por lo que es importante mencionar que si usted quisiera anotar secuencias fragmentadas, samtools podría llegar a fallar al intentar extraer las secuencias de sus genes. Esto se debe a que, al realizar la anotación de los genes con secuencias de proteinas como referencia, para incluir el codón de paro, se busca tres bases después de lo indicado por BLASTx; de modo que si el contig en el que se anotó el gen es más pequeño, no obtendrá la secuencia de este.


# Instalación:

Para clonar el repositorio deberá ejecutar:

```
 git clone https://github.com/marco315176/Rabies_gene_annotator
```

Una vez que haya clonado el repositorio, deberá drigirse al directorio ***bin/*** y otorgarle permiso de ejecución a los scripts:

```
chmod +x *sh
```
También es importante que agregue esta carpeta al PATH en su ***~/.bashrc***:
```
nano ~/.bashrc

export PATH="/$HOME/PATH_TO/Rabies_gene_annotator/bin:$PATH"

source ~/.bashrc

```

# Dependencias necesarias:

Deberá tener instalado los siguientes programas y, de igual forma, deberá agregar los binarios de estos a su PATH:   
***-> BLAST+ para poder ejecutar BLASTx (https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables)  
-> samtools para poder ejecutar samtools faidx (https://www.htslib.org/) y  
-> seqkit (https://bioinf.shenwei.me/seqkit/download/)***

# Preparar la base de datos

Una vez que tenga instaladas las dependencias necesarias, debe ejecutar el script ***RGA_db_dwl.sh*** del siguiente modo:

```
bash RGA_db_dwl.sh
```

Esto descargará las bases de datos en la carpeta $HOME/db/RGA. Si gusta puede agregar la ruta donde se encuentren estos archivos generados a su ~/.bashrc del siguiente modo:

```
nano ~/.bashrc

export Bx_RABV_RGA_PATH="$HOME/db/RGA"

source ~/.bashrc
```


# Uso y ejecución del pipline

Una vez que tenga todos los requisitos, deberá ejecutar el pipeline de este modo:

```
bash rabies_gene_annotator.sh -f FASTA file directory -o OUTPUT directory -p PATH to BLAST DB 

Options:
Usage: /home/bioinfocenasa/bin/rabies_gene_annotator.sh -f FASTA file PATH -o OUTDIR PATH
 -h print help 
 -f FASTA file directory 
 -o OUTPUT directory 
 -p PATH to BLAST database. If you downloades the database by running the RGA_db_dwl.sh script, the path to your database is: /home/bioinfocenasa/db/RGA 

```

# Archivos de salida

Al final de la ejecución del pipline, en el directorio que usted indicó en la opción ***-o***, encontrará dos carpetas: Nucleotidos y Proteinas, en las cuales encontrará las secuencias de nucleotidos y aminoácidos, respectivamente; así como también encontrará el archivo **Annotation_all.tsv** donde estará la información de la anotación de sus secuencias.



