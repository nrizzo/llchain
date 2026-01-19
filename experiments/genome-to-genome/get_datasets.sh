#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

if ! command -v datasets >/dev/null 2>&1
then
    echo "NCBI datasets could not be found!"
    exit 1
fi

human="GCF_000001405.40_GRCh38.p14_genomic.fna"
chimp="GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
arabthaliana="GCF_000001735.4_TAIR10.1_genomic.fna"
arablyrata="GCF_000004255.2_v.1.0_genomic.fna"

for g in $human $chimp $arabthaliana $arablyrata
do
	if [ ! -f "input/${g}.gz" ]
	then
		datasets download genome accession "${g%_${g#*_*_}}"
		
		unzip -p ncbi_dataset.zip "ncbi_dataset/data/${g%_${g#*_*_}}/$g" \
			| gzip -c > "input/$g.gz"
		rm ncbi_dataset.zip
	fi
done
