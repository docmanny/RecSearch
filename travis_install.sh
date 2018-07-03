#!/bin/bash
# Quick script to install some testing dependencies for travisci

cd $HOME
mkdir bin && export PATH="$HOME/bin/:$PATH"
mkdir db && cd db
mkdir blastdb && export BLASTDB="$HOME/db/blastdb"
mkdir blatdb && export BLATDB="$HOME/db/blatdb"
## BLAT
cd blatdb
#rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/myoLuc2/bigZips/myoLuc2.2bit \
#rsync://hgdownload.cse.ucsc.edu/goldenPath/pteVam1/bigZips/pteVam1.2bit \
#rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit ./
touch hg38.2bit, pteVam1.2bit, myoLuc2.2bit
rsync -avzn rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/\
{faToTwoBit,blat/gfClient,blat/gfServer} $HOME/bin/

## BLASTDB Collection
cd ../blastdb

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
Pteropus_vampyrus_transcript_v2.0.fa.nsq

#sudo apt-get install ncbi-blast+
