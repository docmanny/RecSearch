#!/bin/bash
# Quick script to install some testing dependencies for Travis CI

## Create database folders
cd $HOME && echo -e "Changed directory to $HOME"
mkdir db && cd db && echo -e "Changed directory to $HOME/db" || echo -e "COULD NOT CHANGE DIRECTORY TO $HOME/db!!!"
mkdir blastdb && echo -e "Created dir $PWD/blatdb" && export BLASTDB='$HOME/db/blastdb'
mkdir blatdb && echo -e "Created dir $PWD/blatdb" && export BLATDB='$HOME/db/blatdb'

## Add blat databases
cd ~/db/blatdb echo -e "Changed directory to $PWD" || echo -e "COULD NOT CHANGE DIRECTORY TO $HOME/db/blatdb!!!"
#rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/myoLuc2/bigZips/myoLuc2.2bit \
#rsync://hgdownload.cse.ucsc.edu/goldenPath/pteVam1/bigZips/pteVam1.2bit \
#rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit ./ \
# && echo -e "Downloaded all blat .2bit files!" || echo -e "DOWNLOADING BLAT .2bit FILES FAILED!!!"
touch hg38.2bit pteVam1.2bit myoLuc2.2bit \
&& echo -e "Created dummy BLAT files!" || echo -e "FAILED TO MAKE DUMMY BLAT FILES!!!"

## Add blast databases
cd ~/db/blastdb echo -e "Changed directory to $PWD" || echo -e "COULD NOT CHANGE DIRECTORY TO $HOME/db/blastdb!!!"
# Download hg38, mm13 from NCBI and touch the rest to fake them with mocker package
touch Homo_sapiens_genome_GRCh38p9.fa.nhd \
Homo_sapiens_genome_GRCh38p9.fa.nhi \
Homo_sapiens_genome_GRCh38p9.fa.nhr \
Homo_sapiens_genome_GRCh38p9.fa.nin \
Homo_sapiens_genome_GRCh38p9.fa.nog \
Homo_sapiens_genome_GRCh38p9.fa.nsd \
Homo_sapiens_genome_GRCh38p9.fa.nsi \
Homo_sapiens_genome_GRCh38p9.fa.nsq \
Homo_sapiens_protein_GRCh38p9.fa.phd \
Homo_sapiens_protein_GRCh38p9.fa.phi \
Homo_sapiens_protein_GRCh38p9.fa.phr \
Homo_sapiens_protein_GRCh38p9.fa.pin \
Homo_sapiens_protein_GRCh38p9.fa.pog \
Homo_sapiens_protein_GRCh38p9.fa.psd \
Homo_sapiens_protein_GRCh38p9.fa.psi \
Homo_sapiens_protein_GRCh38p9.fa.psq \
Homo_sapiens_transcript_GRCh38p9.fa.nhd \
Homo_sapiens_transcript_GRCh38p9.fa.nhi \
Homo_sapiens_transcript_GRCh38p9.fa.nhr \
Homo_sapiens_transcript_GRCh38p9.fa.nin \
Homo_sapiens_transcript_GRCh38p9.fa.nog \
Homo_sapiens_transcript_GRCh38p9.fa.nsd \
Homo_sapiens_transcript_GRCh38p9.fa.nsi \
Homo_sapiens_transcript_GRCh38p9.fa.nsq \
Myotis_lucifugus_genome_v2.0.fa.nhd \
Myotis_lucifugus_genome_v2.0.fa.nhi \
Myotis_lucifugus_genome_v2.0.fa.nhr \
Myotis_lucifugus_genome_v2.0.fa.nin \
Myotis_lucifugus_genome_v2.0.fa.nog \
Myotis_lucifugus_genome_v2.0.fa.nsd \
Myotis_lucifugus_genome_v2.0.fa.nsi \
Myotis_lucifugus_genome_v2.0.fa.nsq \
Myotis_lucifugus_protein_v2.0.fa.phd \
Myotis_lucifugus_protein_v2.0.fa.phi \
Myotis_lucifugus_protein_v2.0.fa.phr \
Myotis_lucifugus_protein_v2.0.fa.pin \
Myotis_lucifugus_protein_v2.0.fa.pog \
Myotis_lucifugus_protein_v2.0.fa.psd \
Myotis_lucifugus_protein_v2.0.fa.psi \
Myotis_lucifugus_protein_v2.0.fa.psq \
Myotis_lucifugus_transcript_v2.0.fa.nhd \
Myotis_lucifugus_transcript_v2.0.fa.nhi \
Myotis_lucifugus_transcript_v2.0.fa.nhr \
Myotis_lucifugus_transcript_v2.0.fa.nin \
Myotis_lucifugus_transcript_v2.0.fa.nog \
Myotis_lucifugus_transcript_v2.0.fa.nsd \
Myotis_lucifugus_transcript_v2.0.fa.nsi \
Myotis_lucifugus_transcript_v2.0.fa.nsq \
Pteropus_vampyrus_genome_v2.0.fa.nhd \
Pteropus_vampyrus_genome_v2.0.fa.nhi \
Pteropus_vampyrus_genome_v2.0.fa.nhr \
Pteropus_vampyrus_genome_v2.0.fa.nin \
Pteropus_vampyrus_genome_v2.0.fa.nog \
Pteropus_vampyrus_genome_v2.0.fa.nsd \
Pteropus_vampyrus_genome_v2.0.fa.nsi \
Pteropus_vampyrus_genome_v2.0.fa.nsq \
Pteropus_vampyrus_protein_v2.0.fa.phd \
Pteropus_vampyrus_protein_v2.0.fa.phi \
Pteropus_vampyrus_protein_v2.0.fa.phr \
Pteropus_vampyrus_protein_v2.0.fa.pin \
Pteropus_vampyrus_protein_v2.0.fa.pog \
Pteropus_vampyrus_protein_v2.0.fa.psd \
Pteropus_vampyrus_protein_v2.0.fa.psi \
Pteropus_vampyrus_protein_v2.0.fa.psq \
Pteropus_vampyrus_transcript_v2.0.fa.nhd \
Pteropus_vampyrus_transcript_v2.0.fa.nhi \
Pteropus_vampyrus_transcript_v2.0.fa.nhr \
Pteropus_vampyrus_transcript_v2.0.fa.nin \
Pteropus_vampyrus_transcript_v2.0.fa.nog \
Pteropus_vampyrus_transcript_v2.0.fa.nsd \
Pteropus_vampyrus_transcript_v2.0.fa.nsi \
Pteropus_vampyrus_transcript_v2.0.fa.nsq \
&& echo -e "Created dummy BLAST files!" || echo -e "FAILED TO MAKE DUMMY BLAST FILES!!!"

## Add blat and associated files:
cd $HOME && mkdir Downloads || echo -e "Directory $HOME/Downloads already exists"
cd Downloads && echo -e "Changed directory to $HOME/Downloads" || echo -e "COULD NOT CHANGE DIRECTORY TO DOWNLOADS!!!"
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfClient
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/gfServer
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit

mv blat /usr/bin/ && chmod +x /usr/bin/blat
mv gfClient /usr/bin/ && chmod +x /usr/bin/gfClient
mv gfServer /usr/bin/ && chmod +x /usr/bin/gfServer
mv twoBitToFa /usr/bin/ && chmod +x /usr/bin/twoBitToFa
mv faToTwoBit /usr/bin/ && chmod +x /usr/bin/faToTwoBit

## Download RefSeq transcript files and use it to make .2bit files:
cd $HOME && mkdir seq
curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/GCF_000001405.37_GRCh38.p11_rna.fna.gz' -o ./HomSapGRCH38_transcript.fa.gz
gunzip -d HomSapGRCH38_transcript.fa.gz && rm HomSapGRCH38_transcript.fa.gz
faToTwoBit HomSapGRCH38_transcript.fa Homo_sapiens_transcript.2bit && mv Homo_sapiens_transcript.2bit ~/db/blatdb/

