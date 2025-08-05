# Rabies_gene_annotator
Pipline escrito en bash para la anotación y obtención de los 5 genes del virus de la rabia

Importante:
Tenga en cuenta que sus secuencias deben estar nombradas algo así; de lo contrario el script fallará al intentar reconocer y formular los ID's: RABV-2025-Rab90.fa

Si lo desea puede hacer los cambios pertinentes en las lineas de cada ciclo "for", para señalar cuál es el ID de las secuencias: 
ID=$(basename ${assembly} | cut -d '-' -f '3' | cut -d '.' -f '1')
inf=$(basename ${assembly} | cut -d '-' -f '2')
newID=${ID}-${inf}
