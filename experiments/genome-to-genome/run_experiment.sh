#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

if ! command -v seqtk >/dev/null 2>&1
then
        echo "ERROR: seqtk not found!"; exit 1
fi

mkdir output
cd output

clcviz=$thisfolder/../../clc-viz
usrbintime="/usr/bin/time -f"
usrbintimeoptions="%e total time elapsed (s)\n%M maxresident k"

echo "# human vs chimp (chr1 only)"
mkfifo human.fa.gz
refs=$thisfolder/output/human.fa.gz
mkfifo chimp.fa.gz
queries=$thisfolder/output/chimp.fa.gz
clcvizmode="--mode global"
clcvizseed="-a MUM -l 20"

# run ChainX* (via clc-viz)
seqtk subseq $thisfolder/input/GCF_000001405.40_GRCh38.p14_genomic.fna.gz <(echo NC_000001.11) > human.fa.gz &
seqtk subseq $thisfolder/input/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz <(echo NC_072398.2) > chimp.fa.gz &
$usrbintime "$usrbintimeoptions" $clcviz --chainx $clcvizmode $clcvizseed -t $refs -q $queries \
        >> human_mum_chainx 2>> human_mum_chainx
# run ChainX*-opt (via clc-viz)
seqtk subseq $thisfolder/input/GCF_000001405.40_GRCh38.p14_genomic.fna.gz <(echo NC_000001.11) > human.fa.gz &
seqtk subseq $thisfolder/input/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz <(echo NC_072398.2) > chimp.fa.gz &
$usrbintime "$usrbintimeoptions" $clcviz --chainx-opt $clcvizmode $clcvizseed -t $refs -q $queries \
        >> human_mum_chainx-opt 2>> human_mum_chainx-opt
# run clc-viz
seqtk subseq $thisfolder/input/GCF_000001405.40_GRCh38.p14_genomic.fna.gz <(echo NC_000001.11) > human.fa.gz &
seqtk subseq $thisfolder/input/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz <(echo NC_072398.2) > chimp.fa.gz &
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -t $refs -q $queries \
        >> human_mum_clc-viz 2>> human_mum_clc-viz

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "anchored edit distance" human_mum_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" human_mum_clc-viz    | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_human_mum_headers
echo "$clcvizseed" >> stats_human_mum_headers
for t in "human_mum_chainx" "human_mum_chainx-opt" "human_mum_clc-viz"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_human_mum_headers \
	stats_time_human_mum_chainx     stats_space_human_mum_chainx \
	stats_time_human_mum_chainx-opt stats_space_human_mum_chainx-opt \
	stats_time_human_mum_clc-viz    stats_space_human_mum_clc-viz \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$clc-viz") - | column -t -s'$'

# time per input anchor
outputpng="times_human_mum.png"
echo -n "Plotting individual times in output/$outputpng ..."
for t in "human_mum_chainx" "human_mum_chainx-opt" "human_mum_clc-viz"
do
	# anchor, total time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f28) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_human_mum_chainx', 'times_human_mum_chainx-opt', 'times_human_mum_clc-viz'"
echo "done."

echo "# Arabidopsis thaliana vs Arabidopsis lyrata"
refs=$thisfolder/input/GCF_000001735.4_TAIR10.1_genomic.fna.gz
queries=$thisfolder/input/GCF_000004255.2_v.1.0_genomic.fna.gz
clcvizmode="--mode global"
clcvizseed="-a MUM -l 20"

# run ChainX*
$usrbintime "$usrbintimeoptions" $clcviz --chainx $clcvizmode $clcvizseed -t $refs -q $queries \
	>> plant_mum_chainx 2>> plant_mum_chainx
# run ChainX*-opt
$usrbintime "$usrbintimeoptions" $clcviz --chainx-opt $clcvizmode $clcvizseed -t $refs -q $queries \
	>> plant_mum_chainx-opt 2>> plant_mum_chainx-opt
# run clc-viz
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -t $refs -q $queries \
	>> plant_mum_clc-viz 2>> plant_mum_clc-viz

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "anchored edit distance" plant_mum_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" plant_mum_clc-viz    | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_plant_mum_headers
echo "$clcvizseed" >> stats_plant_mum_headers
for t in "plant_mum_chainx" "plant_mum_chainx-opt" "plant_mum_clc-viz"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_plant_mum_headers \
	stats_time_plant_mum_chainx     stats_space_plant_mum_chainx \
	stats_time_plant_mum_chainx-opt stats_space_plant_mum_chainx-opt \
	stats_time_plant_mum_clc-viz    stats_space_plant_mum_clc-viz \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$clc-viz") - | column -t -s'$'

# time per input anchor
outputpng="times_plant_mum.png"
for t in "plant_mum_clc-viz" "plant_mum_chainx" "plant_mum_chainx-opt"
do
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f28) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_plant_mum_chainx', 'times_plant_mum_chainx-opt', 'times_plant_mum_clc-viz'"
echo "done."

echo "# HPRC human vs HPRC human (chr19 only)"
# https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/annotation/chrom_assignment/chrom_alias_hprc_r2_v1.0.index.csv
( \
	seqtk subseq $thisfolder/input/HG00097_hap1_hprc_r2_v1.0.1.fa.gz <(echo "HG00097#1#CM094070.1") ;
	seqtk subseq $thisfolder/input/HG00097_hap2_hprc_r2_v1.0.1.fa.gz <(echo "HG00097#2#CM094085.1") ;
	seqtk subseq $thisfolder/input/HG00272_hap1_hprc_r2_v1.0.1.fa.gz <(echo "HG00272#1#CM094199.1") ;
	seqtk subseq $thisfolder/input/HG00272_hap2_hprc_r2_v1.0.1.fa.gz <(echo "HG00272#2#CM094225.1") ;
	seqtk subseq $thisfolder/input/HG00280_hap1_hprc_r2_v1.0.1.fa.gz <(echo "HG00280#1#CM087157.1") ;
	seqtk subseq $thisfolder/input/HG00280_hap2_hprc_r2_v1.0.1.fa.gz <(echo "HG00280#2#CM087169.1") ;
	seqtk subseq $thisfolder/input/HG00408_pat_hprc_r2_v1.0.1.fa.gz  <(echo "HG00408#1#CM085963.1") ;
	seqtk subseq $thisfolder/input/HG00408_mat_hprc_r2_v1.0.1.fa.gz  <(echo "HG00408#2#CM085976.1") ;
) | gzip -c > hprc_chr19.fa.gz

queries=hprc_chr19.fa.gz
clcvizmode="--mode global --all-to-all"
clcvizseed="-a MUM -l 20"

# run ChainX*
$usrbintime "$usrbintimeoptions" $clcviz --chainx $clcvizmode $clcvizseed -q $queries \
	>> hprc_mum_chainx.phylip 2>> hprc_mum_chainx
# run ChainX*-opt
$usrbintime "$usrbintimeoptions" $clcviz --chainx-opt $clcvizmode $clcvizseed -q $queries \
	>> hprc_mum_chainx-opt.phylip 2>> hprc_mum_chainx-opt
# run clc-viz
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -q $queries \
	>> hprc_mum_clc-viz.phylip 2>> hprc_mum_clc-viz

echo -n "Checking if the optimal chaining cost (ChainX-opt vs clc-viz) differs..."
check=$(diff \
	<(grep "anchored edit distance" hprc_mum_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" hprc_mum_clc-viz    | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_hprc_mum_headers
echo "$clcvizseed" >> stats_hprc_mum_headers
for t in "hprc_mum_chainx" "hprc_mum_chainx-opt" "hprc_mum_clc-viz"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_hprc_mum_headers \
	stats_time_hprc_mum_chainx     stats_space_hprc_mum_chainx \
	stats_time_hprc_mum_chainx-opt stats_space_hprc_mum_chainx-opt \
	stats_time_hprc_mum_clc-viz    stats_space_hprc_mum_clc-viz \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$clc-viz") - | column -t -s'$'

# time per input anchor
outputpng="times_hprc_mum.png"
echo -n "Plotting individual times in output/$outputpng ..."
for t in "hprc_mum_chainx" "hprc_mum_chainx-opt" "hprc_mum_clc-viz"
do
	# anchor, total time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f26) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_hprc_mum_chainx', 'times_hprc_mum_chainx-opt', 'times_hprc_mum_clc-viz'"
echo "done."
