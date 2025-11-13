#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder
mkdir output
cd output

clcviz=$thisfolder/../../clc-viz
chainx=$thisfolder/../../ext/ChainX/chainX
usrbintime="/usr/bin/time -f"
usrbintimeoptions="%e total time elapsed (s)\n%M maxresident k"

echo "# MUMs (l >= 20)"
ref=$thisfolder/../../ext/ChainX/data/human/chm13v2.0_concat.fa.gz
queries=$thisfolder/../../ext/ChainX/data/human/sample_100k.fa.gz
clcvizmode="--mode semiglobal"
chainxmode="-m sg"
clcvizseed="-a MUM -l 20"
chainxseed="-a MUM -l 20"

# run ChainX*
$usrbintime "$usrbintimeoptions" $chainx $chainxmode $chainxseed -t $ref -q $queries \
	>> human_mum_chainx 2>> human_mum_chainx
# run ChainX*-opt
$usrbintime "$usrbintimeoptions" $chainx $chainxmode $chainxseed -t $ref -q $queries --optimal \
	>> human_mum_chainx-opt 2>> human_mum_chainx-opt
# run clc-viz
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -t $ref -q $queries \
	>> human_mum_clc-viz 2>> human_mum_clc-viz

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "distance =" human_mum_chainx-opt | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "anchored edit distance" human_mum_clc-viz | cut -d' ' -f16) ; exit 0)
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
for t in "human_mum_chainx" "human_mum_chainx-opt"
do
	paste \
		<(grep "count of anchors" $t | cut -d' ' -f9 | tr -d ",") \
		<(grep "distance computation finished" $t | cut -d' ' -f6 | tr -d "(") \
		> times_$t
done
for t in "human_mum_clc-viz"
do
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f26) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_human_mum_chainx', 'times_human_mum_chainx-opt', 'times_human_mum_clc-viz'"
echo "done."

echo "# MEMs (l >= 75)"
clcvizseed="-a MEM -l 75"
chainxseed="-a MEM -l 75"

# run ChainX*
$usrbintime "$usrbintimeoptions" $chainx $chainxmode $chainxseed -t $ref -q $queries \
	>> human_mem_chainx 2>> human_mem_chainx
# run ChainX*-opt
$usrbintime "$usrbintimeoptions" $chainx $chainxmode $chainxseed -t $ref -q $queries --optimal \
	>> human_mem_chainx-opt 2>> human_mem_chainx-opt
# run clc-viz
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -t $ref -q $queries \
	>> human_mem_clc-viz 2>> human_mem_clc-viz

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "distance =" human_mem_chainx-opt | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "anchored edit distance" human_mem_clc-viz | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_human_mem_headers
echo "$clcvizseed" >> stats_human_mem_headers
for t in "human_mem_chainx" "human_mem_chainx-opt" "human_mem_clc-viz"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_human_mem_headers \
	stats_time_human_mem_chainx     stats_space_human_mem_chainx \
	stats_time_human_mem_chainx-opt stats_space_human_mem_chainx-opt \
	stats_time_human_mem_clc-viz    stats_space_human_mem_clc-viz \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$clc-viz") - | column -t -s'$'

# time per input anchor
outputpng="times_human_mem.png"
echo -n "Plotting individual times in output/$outputpng ..."
for t in "human_mem_chainx" "human_mem_chainx-opt"
do
	paste \
		<(grep "count of anchors" $t | cut -d' ' -f9 | tr -d ",") \
		<(grep "distance computation finished" $t | cut -d' ' -f6 | tr -d "(") \
		> times_$t
done
for t in "human_mem_clc-viz"
do
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f26) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_human_mem_chainx', 'times_human_mem_chainx-opt', 'times_human_mem_clc-viz'"
echo "done."
