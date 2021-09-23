# download the miRNA  binding sits in 3UTR from the miRwalk database
wget http://mirwalk.umm.uni-heidelberg.de/download/hsa_miRWalk_3UTR.zip
gzip -d hsa_miRWalk_3UTR.zip

#download the RBP binding sits in 3UTR from beRBP database
wget http://bioinfo.vanderbilt.edu/beRBP/download/beRBP-G_targets_26RBPs.7z
7z x beRBP-G_targets_26RBPs.7z
grep -H ">" *.fa > all_RBP_with_RBPnames.txt

# download the UCSC-ensemble id files from UCSC( the version of gene id in the beRBP files is hg19)
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownToEnsembl.txt.gz
gzip -d knownToEnsembl.txt.gz

# download cancer-gene relationship from cancermine datavase
wget https://zenodo.org/record/5363661/files/cancermine_collated.tsv?download=1

# download protein-protein interaction from STRING database
wget https://stringdb-static.org/download/protein.physical.links.full.v11.5/9606.protein.physical.links.full.v11.5.txt.gz

# download m6A data from the REPIC database
wget https://repicmod.uchicago.edu/repic/data/download/m6A=sites=species=human=hg38.txt.gz
gzip -d m6A=sites=species=human=hg38.txt.gz
grep utr3 m6A=sites=species=human=hg38.txt > utr3m6a.txt
