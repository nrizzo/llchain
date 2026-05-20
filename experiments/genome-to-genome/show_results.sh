#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder
cd output

echo "# human vs chimp (chr1 only)"
llchainmode="--mode global"
llchainseed="-a MUM -l 20"

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "anchored edit distance" human_mum_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" human_mum_llchain    | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_human_mum_headers
echo "$llchainseed" >> stats_human_mum_headers
for t in "human_mum_chainx" "human_mum_chainx-opt" "human_mum_llchain"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_human_mum_headers \
	stats_time_human_mum_chainx     stats_space_human_mum_chainx \
	stats_time_human_mum_chainx-opt stats_space_human_mum_chainx-opt \
	stats_time_human_mum_llchain    stats_space_human_mum_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

# time per input anchor
outputpng="times_human_mum.png"
echo -n "Plotting individual times in output/$outputpng ..."
for t in "human_mum_chainx" "human_mum_chainx-opt" "human_mum_llchain"
do
	# anchor, total time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f28) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_human_mum_chainx', 'times_human_mum_chainx-opt', 'times_human_mum_llchain'"
echo "done."

echo "# Arabidopsis thaliana vs Arabidopsis lyrata"
llchainseed="-a MUM -l 20"

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "anchored edit distance" plant_mum_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" plant_mum_llchain    | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_plant_mum_headers
echo "$llchainseed" >> stats_plant_mum_headers
for t in "plant_mum_chainx" "plant_mum_chainx-opt" "plant_mum_llchain"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_plant_mum_headers \
	stats_time_plant_mum_chainx     stats_space_plant_mum_chainx \
	stats_time_plant_mum_chainx-opt stats_space_plant_mum_chainx-opt \
	stats_time_plant_mum_llchain    stats_space_plant_mum_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

# time per input anchor
outputpng="times_plant_mum.png"
for t in "plant_mum_llchain" "plant_mum_chainx" "plant_mum_chainx-opt"
do
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f28) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_plant_mum_chainx', 'times_plant_mum_chainx-opt', 'times_plant_mum_llchain'"
echo "done."

echo "# HPRC human vs HPRC human (chr19 only)"
llchainmode="--mode global --all-to-all"
llchainseed="-a MUM -l 20"

echo -n "Checking if the optimal chaining cost (ChainX-opt vs llchain) differs..."
check=$(diff \
	<(grep "anchored edit distance" hprc_mum_chainx-opt | cut -d' ' -f16) \
	<(grep "anchored edit distance" hprc_mum_llchain    | cut -d' ' -f16) ; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

# total time, memory
echo "seeds" > stats_hprc_mum_headers
echo "$llchainseed" >> stats_hprc_mum_headers
for t in "hprc_mum_chainx" "hprc_mum_chainx-opt" "hprc_mum_llchain"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t
done
paste -d'$' stats_hprc_mum_headers \
	stats_time_hprc_mum_chainx     stats_space_hprc_mum_chainx \
	stats_time_hprc_mum_chainx-opt stats_space_hprc_mum_chainx-opt \
	stats_time_hprc_mum_llchain    stats_space_hprc_mum_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

# time per input anchor
outputpng="times_hprc_mum.png"
echo -n "Plotting individual times in output/$outputpng ..."
for t in "hprc_mum_chainx" "hprc_mum_chainx-opt" "hprc_mum_llchain"
do
	# anchor, total time
	paste \
		<(grep "querying" $t | cut -d' ' -f14) \
		<(grep "querying" $t | cut -d' ' -f26) \
		> times_$t
done
gnuplot -persist -e "set term pngcairo; set xlabel \"anchors\"; set ylabel \"seconds\"; set title \"Seeding+chaining time\"; set output '$outputpng'; plot 'times_hprc_mum_chainx', 'times_hprc_mum_chainx-opt', 'times_hprc_mum_llchain'"
echo "done."
