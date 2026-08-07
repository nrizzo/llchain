#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder
cd output

echo "# MUMs (l >= 20)"
llchainseed="-a MUM -l 20"

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep -v "^\[llchain\]" human_mum_chainx-opt | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f3) \
	<(grep -v "^\[llchain\]"    human_mum_llchain | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f3) \
	; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

echo "seeds" > stats_human_mum_headers
echo "$llchainseed" >> stats_human_mum_headers
for t in "human_mum_chainx" "human_mum_chainx-opt" "human_mum_llchain"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	echo "avg chaining time (s)" > stats_avg_chaining_time_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t

	# anchor, total time
	grep -v "^\[llchain\]" $t | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f4,10 \
		> total_times_$t

	# anchor, seeding time
	grep -v "^\[llchain\]" $t | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f4,6 \
		> seeding_times_$t

	# anchor, preprocess + chain + postprocess time
	grep -v "^\[llchain\]" $t | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f4,7,8,9 \
		| awk 'BEGIN {OFS="\t"} {print $1,$2+$3+$4}' \
		> chaining_times_$t

	# avg chaining time
	awk 'BEGIN {total = 0; n = 0} {total += $2; n += 1} END {print total / n}' chaining_times_$t \
		>> stats_avg_chaining_time_$t
done
paste -d'$' stats_human_mum_headers \
	stats_time_human_mum_chainx     stats_avg_chaining_time_human_mum_chainx     stats_space_human_mum_chainx \
	stats_time_human_mum_chainx-opt stats_avg_chaining_time_human_mum_chainx-opt stats_space_human_mum_chainx-opt \
	stats_time_human_mum_llchain    stats_avg_chaining_time_human_mum_llchain    stats_space_human_mum_llchain \
	| cat <(echo -e "\$ChainX*\$\$\$ChainX-opt*\$\$\$llchain") - | column -t -s'$'

echo "# MEMs (l >= 50)"
llchainseed="-a MEM -l 50"

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep -v "^\[llchain\]" human_mem_chainx-opt | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f3) \
	<(grep -v "^\[llchain\]"    human_mem_llchain | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f3) \
	; exit 0)
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

echo "seeds" > stats_human_mem_headers
echo "$llchainseed" >> stats_human_mem_headers
for t in "human_mem_chainx" "human_mem_chainx-opt" "human_mem_llchain"
do
	echo "time (s)" > stats_time_$t
	echo "space (kb)" > stats_space_$t
	echo "avg chaining time (s)" > stats_avg_chaining_time_$t
	grep "total time elapsed" $t | cut -d' ' -f1 >> stats_time_$t
	grep "maxresident" ${t} | cut -d' ' -f1 >> stats_space_$t

	# anchor, total time
	grep -v "^\[llchain\]" $t | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f4,10 \
		> total_times_$t

	# anchor, seeding time
	grep -v "^\[llchain\]" $t | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f4,6 \
		> seeding_times_$t

	# anchor, preprocess + chain + postprocess time
	grep -v "^\[llchain\]" $t | grep -v "^#" | grep -v "maxresident" | grep -v "elapsed" | cut -f4,7,8,9 \
		| awk 'BEGIN {OFS="\t"} {if ($1 > 0) print $1,$2+$3+$4}' \
		> chaining_times_$t

	# avg chaining time
	awk 'BEGIN {total = 0; n = 0} {total += $2; n += 1} END {print total / n}' chaining_times_$t \
		>> stats_avg_chaining_time_$t
done
paste -d'$' stats_human_mem_headers \
	stats_time_human_mem_chainx     stats_avg_chaining_time_human_mem_chainx     stats_space_human_mem_chainx \
	stats_time_human_mem_chainx-opt stats_avg_chaining_time_human_mem_chainx-opt stats_space_human_mem_chainx-opt \
	stats_time_human_mem_llchain    stats_avg_chaining_time_human_mem_llchain    stats_space_human_mem_llchain \
	| cat <(echo -e "\$ChainX*\$\$\$ChainX-opt*\$\$\$llchain") - | column -t -s'$'

echo "Computing output/plot.svg..."
python3 $thisfolder/plot_times_log_log.py
