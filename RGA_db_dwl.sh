#!/bin/bash

#---------------------------
# Database download for RGA
#---------------------------

echo -e "\033[42m================================================\033[m"
echo -e "\033[42m========== Starting database download ==========\033[m"
echo -e "\033[42m================================================\033[m""\n"

echo -e "\033[4;33m========== Creating database in $HOME/db/RGA ==========\033[0m\n"

mkdir -p $HOME/db/RGA
cd $HOME/db/RGA

#--------
# N prot
#--------

esearch -db protein -query "Lyssavirus rabies [organism] AND nucleoprotein [Title]" \
       | efetch -format fasta > Nucleoprot_RABV.faa

seqkit seq -g -m 450 -M 450 Nucleoprot_RABV.faa > Nucleoprot_RABV.fa

echo -e "\033[0;36m========== N gene database created ==========\033[0m\n"

#--------
# P prot
#--------

esearch -db protein -query "Lyssavirus rabies [organism] AND phosphoprotein [Title]" \
       | efetch -format fasta > Phosphoprot_RABV.faa

seqkit seq -g -m 297 -M 297 Phosphoprot_RABV.faa > Phosphoprot_RABV.fa

echo -e "\033[0;36m========== P gene database created ==========\033[0m\n"

#----------
# Mtx prot
#----------

esearch -db protein -query "Lyssavirus rabies [organism] AND matrix protein [Title]" \
       | efetch -format fasta > Mtxprot_RABV.faa

seqkit seq -g -m 202 -M 202 Mtxprot_RABV.faa > Mtxprot_RABV.fa

echo -e "\033[0;36m========== M gene database created ==========\033[0m\n"

#---------
# G prot
#---------

esearch -db protein -query "Lyssavirus rabies [organism] AND glycoprotein [Title]" \
       | efetch -format fasta > Glicoprot_RABV.faa

seqkit seq -g -m 524 -M 524 Glicoprot_RABV.faa > Glicoprot_RABV.fa

echo -e "\033[0;36m========== G gene database created ==========\033[0m\n"

#--------
# L prot
#--------

esearch -db protein -query "Lyssavirus rabies [organism] AND polymerase [Title]" \
       | efetch -format fasta > ProtL_RABV.faa

seqkit seq -g -m 2127 -M 2127 ProtL_RABV.faa > ProtL_RABV.fa

echo -e "\033[0;36m========== L gene database created ==========\033[0m\n"

chmod +x ./*.fa
rm ./*.faa


#--------------------------------
# Creating a database with BLAST
#--------------------------------

makeblastdb -in Nucleoprot_RABV.fa -dbtype prot -out ./Nprot_RABV_db #For nucleoprotein (N)
makeblastdb -in Phosphoprot_RABV.fa -dbtype prot -out ./Pprot_RABV_db #For phosphoprotein (P)
makeblastdb -in Mtxprot_RABV.fa -dbtype prot -out ./Mtxprot_RABV_db #For matrix protein (M)
makeblastdb -in Glicoprot_RABV.fa -dbtype prot -out ./Gprot_RABV_db #For glycoprotein (G)
makeblastdb -in ProtL_RABV.fa -dbtype prot -out ./Lprot_RABV_db #For polymerase protein (L)

 echo -e "\033[46m========== databases created ==========\033[m"
