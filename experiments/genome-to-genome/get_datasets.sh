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

hprc_assemblies=(
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00097/assemblies/freeze_2/HG00097_hap2_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00097/assemblies/freeze_2/HG00097_hap1_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00272/assemblies/freeze_2/HG00272_hap2_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00272/assemblies/freeze_2/HG00272_hap1_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00280/assemblies/freeze_2/HG00280_hap2_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00280/assemblies/freeze_2/HG00280_hap1_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00408/assemblies/freeze_2/HG00408_mat_hprc_r2_v1.0.1.fa.gz"
	"https://human-pangenomics.s3.amazonaws.com/submissions/DC27718F-5F38-43B0-9A78-270F395F13E8--INT_ASM_PRODUCTION/HG00408/assemblies/freeze_2/HG00408_pat_hprc_r2_v1.0.1.fa.gz"
	)

for link in ${hprc_assemblies[@]}
do
	if [ ! -f "input/${link##*/}" ]
	then
		wget "$link" -O "input/${link##*/}"
	fi
done
