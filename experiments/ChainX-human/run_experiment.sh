#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder
mkdir output
cd output

clcviz=$thisfolder/../../clc-viz
usrbintime="/usr/bin/time -f"
usrbintimeoptions="%e total time elapsed (s)\n%M maxresident k"

echo "# MUMs (l >= 20)"
ref=$thisfolder/../../ext/ChainX/data/human/chm13v2.0_concat.fa.gz
queries=$thisfolder/../../ext/ChainX/data/human/sample_100k.fa.gz
clcvizmode="--mode semiglobal"
clcvizseed="-a MUM -l 20"

# run ChainX* (via clc-viz)
$usrbintime "$usrbintimeoptions" $clcviz --chainx $clcvizmode $clcvizseed -t $ref -q $queries \
        >> human_mum_chainx 2>> human_mum_chainx
# run ChainX*-opt (via clc-viz)
$usrbintime "$usrbintimeoptions" $clcviz --chainx-opt $clcvizmode $clcvizseed -t $ref -q $queries \
        >> human_mum_chainx-opt 2>> human_mum_chainx-opt
# run clc-viz
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -t $ref -q $queries \
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
echo -n "Plotting the times per anchor..."
for t in "human_mum_chainx" "human_mum_chainx-opt" "human_mum_clc-viz"
do
	# anchor, total time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f28) \
		> total_times_$t

	# anchor, seeding time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f20) \
		> seeding_times_$t

	# anchor, preprocess + chain + postprocess time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f22) \
		<(grep "querying" $t | cut -d' ' -f24) \
		<(grep "querying" $t | cut -d' ' -f26) \
		| awk 'BEGIN {OFS="\t"} {print $1,$2+$3+$4}' \
		> chaining_times_$t
done

outputpng="total_times_human_mum.png"
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Total time\"; set output '$outputpng'; plot 'total_times_human_mum_chainx', 'total_times_human_mum_chainx-opt', 'total_times_human_mum_clc-viz'"

outputpng="seeding_times_human_mum.png"
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding time\"; set output '$outputpng'; plot 'seeding_times_human_mum_chainx', 'seeding_times_human_mum_chainx-opt', 'seeding_times_human_mum_clc-viz'"

outputpng="chaining_times_human_mum.png"
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Preprocess+chaining+postprocess time\"; set output '$outputpng'; plot 'chaining_times_human_mum_chainx', 'chaining_times_human_mum_chainx-opt', 'chaining_times_human_mum_clc-viz'"
echo "done."

echo "# MEMs (l >= 50)"
clcvizseed="-a MEM -l 50"
chainxseed="-a MEM -l 50"

# run ChainX*
$usrbintime "$usrbintimeoptions" $clcviz --chainx $clcvizmode $clcvizseed -t $ref -q $queries \
	>> human_mem_chainx 2>> human_mem_chainx
# run ChainX*-opt
$usrbintime "$usrbintimeoptions" $clcviz --chainx-opt $clcvizmode $clcvizseed -t $ref -q $queries \
	>> human_mem_chainx-opt 2>> human_mem_chainx-opt
# run clc-viz
$usrbintime "$usrbintimeoptions" $clcviz $clcvizmode $clcvizseed -t $ref -q $queries \
	>> human_mem_clc-viz 2>> human_mem_clc-viz

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "anchored edit distance" human_mem_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" human_mem_clc-viz    | cut -d' ' -f16) ; exit 0)
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

echo -n "Plotting the times per anchor..."
for t in "human_mem_chainx" "human_mem_chainx-opt" "human_mem_clc-viz"
do
	# anchor, total time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f28) \
		> total_times_$t

	# anchor, seeding time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f20) \
		> seeding_times_$t

	# anchor, preprocess + chain + postprocess time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f22) \
		<(grep "querying" $t | cut -d' ' -f24) \
		<(grep "querying" $t | cut -d' ' -f26) \
		| awk 'BEGIN {OFS="\t"} {print $1,$2+$3+$4}' \
		> chaining_times_$t
done

outputpng="total_times_human_mem.png"
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Total time\"; set output '$outputpng'; plot 'total_times_human_mem_chainx', 'total_times_human_mem_chainx-opt', 'total_times_human_mem_clc-viz'"

outputpng="seeding_times_human_mem.png"
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding time\"; set output '$outputpng'; plot 'seeding_times_human_mem_chainx', 'seeding_times_human_mem_chainx-opt', 'seeding_times_human_mem_clc-viz'"

outputpng="chaining_times_human_mem.png"
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Preprocess+chaining+postprocess time\"; set output '$outputpng'; plot 'chaining_times_human_mem_chainx', 'chaining_times_human_mem_chainx-opt', 'chaining_times_human_mem_clc-viz'"
echo "done."
