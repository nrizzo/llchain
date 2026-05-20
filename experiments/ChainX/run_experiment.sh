#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder

llchain=$thisfolder/../../llchain
chainx=$thisfolder/../../ext/ChainX/chainX
usrbintime="/usr/bin/time -f"
usrbintimeoptions="%M maxresident k"

mkdir output && cd output

echo "# Global sequence comparison"
ref="$thisfolder/../../ext/ChainX/data/time_global/Chromosome_2890043_3890042_0.fasta"
llchainmode="--mode global"
chainxmode="-m g"

# run ChainX*
for m in 70 80 90 94 97 99
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/time_global/mutated_${m}_perc.fasta" \
		>> time_global_chainx 2>> time_global_chainx
done

# run ChainX-opt*
for m in 70 80 90 94 97 99
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/time_global/mutated_${m}_perc.fasta" \
		--optimal >> time_global_chainx_opt 2>> time_global_chainx_opt
done

# run llchain
for m in 70 80 90 94 97 99
do
	$usrbintime "$usrbintimeoptions" $llchain $llchainmode -t $ref -q "$thisfolder/../../ext/ChainX/data/time_global/mutated_${m}_perc.fasta" \
		>> time_global_llchain 2>> time_global_llchain
done

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "distance =" time_global_chainx     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "distance =" time_global_chainx_opt | cut -d'=' -f2 | cut -d')' -f2 | tr -d " "))
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
check=$(diff \
	<(grep "distance =" time_global_chainx     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "anchored edit distance" time_global_llchain | cut -d' ' -f16))
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

echo "MUMs" > stats_time_global_mums
grep "count of anchors" time_global_chainx | cut -d' ' -f9 | tr -d "," | tac >> stats_time_global_mums

# time, memory, iterations
for t in time_global_chainx time_global_chainx_opt
do
	echo "time (s)" > stats_time_$t
	grep "distance computation finished" $t | cut -d'(' -f2 | cut -d' ' -f1 | tac >> stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "maxresident" $t | cut -d' ' -f1 | tac >> stats_space_$t
done
echo "time (s)" > stats_time_time_global_llchain
grep "anchored edit distance" time_global_llchain | cut -d' ' -f28 | tr -d "s" | tac >> stats_time_time_global_llchain
echo "space (kb)" > stats_space_time_global_llchain
grep "maxresident" time_global_llchain | cut -d' ' -f1 | tac >> stats_space_time_global_llchain

paste -d'$' stats_time_global_mums \
	stats_time_time_global_chainx             stats_space_time_global_chainx \
	stats_time_time_global_chainx_opt         stats_space_time_global_chainx_opt  \
	stats_time_time_global_llchain            stats_space_time_global_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

echo "# Semiglobal sequence comparison"
ref="$thisfolder/../../ext/ChainX/data/time_semiglobal/e_coli_DH1_reference.fasta"
llchainmode="--mode semiglobal"
chainxmode="-m sg"

# run ChainX*
for m in 70 80 90 94 97
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/time_semiglobal/mutated_${m}_perc.fasta" \
		>> time_semiglobal_chainx 2>> time_semiglobal_chainx
done

# run ChainX-opt*
for m in 70 80 90 94 97
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/time_semiglobal/mutated_${m}_perc.fasta" \
		--optimal >> time_semiglobal_chainx_opt 2>> time_semiglobal_chainx_opt
done

# run llchain
for m in 70 80 90 94 97
do
	$usrbintime "$usrbintimeoptions" $llchain $llchainmode -t $ref -q "$thisfolder/../../ext/ChainX/data/time_semiglobal/mutated_${m}_perc.fasta" \
		>> time_semiglobal_llchain 2>> time_semiglobal_llchain
done

echo -n "Checking if the optimal chaining cost differs..."
check=$(diff \
	<(grep "distance =" time_semiglobal_chainx     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "distance =" time_semiglobal_chainx_opt | cut -d'=' -f2 | cut -d')' -f2 | tr -d " "))
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
check=$(diff \
	<(grep "distance =" time_semiglobal_chainx     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
	<(grep "anchored edit distance" time_semiglobal_llchain | cut -d' ' -f16))
if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
echo " done (no difference)."

echo "MUMs" > stats_time_semiglobal_mums
grep "count of anchors" time_semiglobal_chainx | cut -d' ' -f9 | tr -d "," | tac >> stats_time_semiglobal_mums

# time, memory, iterations
for t in time_semiglobal_chainx time_semiglobal_chainx_opt
do
	echo "time (s)" > stats_time_$t
	grep "distance computation finished" $t | cut -d'(' -f2 | cut -d' ' -f1 | tac >> stats_time_$t
	echo "space (kb)" > stats_space_$t
	grep "maxresident" $t | cut -d' ' -f1 | tac >> stats_space_$t
done
echo "time (s)" > stats_time_time_semiglobal_llchain
grep "anchored edit distance" time_semiglobal_llchain | cut -d' ' -f28 | tr -d "s" | tac >> stats_time_time_semiglobal_llchain
echo "space (kb)" > stats_space_time_semiglobal_llchain
grep "maxresident" time_semiglobal_llchain | cut -d' ' -f1 | tac >> stats_space_time_semiglobal_llchain

paste -d'$' stats_time_semiglobal_mums \
	stats_time_time_semiglobal_chainx             stats_space_time_semiglobal_chainx \
	stats_time_time_semiglobal_chainx_opt         stats_space_time_semiglobal_chainx_opt  \
	stats_time_time_semiglobal_llchain            stats_space_time_semiglobal_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

echo "# Global sequence comparison, 100 queries"
usrbintimeoptions="%e total time elapsed (s)\n%M maxresident k"
ref="$thisfolder/../../ext/ChainX/data/correlation_global/Chromosome_2890043_3890042_0.fasta"
llchainmode="--mode global"
chainxmode="-m g"

# run ChainX*
for m in "75_80" "80_90" "90_100"
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/correlation_global/count100_mutated_$m.fa.gz" \
		>> correlation_global_chainx_$m 2>> correlation_global_chainx_$m
done

# run ChainX-opt*
for m in "75_80" "80_90" "90_100"
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/correlation_global/count100_mutated_$m.fa.gz" \
		--optimal >> correlation_global_chainx_opt_$m 2>> correlation_global_chainx_opt_$m
done

# run llchain
for m in "75_80" "80_90" "90_100"
do
	$usrbintime "$usrbintimeoptions" $llchain $llchainmode -t $ref -q "$thisfolder/../../ext/ChainX/data/correlation_global/count100_mutated_$m.fa.gz" \
		>> correlation_global_llchain_$m 2>> correlation_global_llchain_$m
done

echo -n "Checking if the optimal chaining cost differs..."
for m in "75_80" "80_90" "90_100"
do
	check=$(diff \
		<(grep "distance =" correlation_global_chainx_$m     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
		<(grep "distance =" correlation_global_chainx_opt_$m | cut -d'=' -f2 | cut -d')' -f2 | tr -d " "))
	if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
	check=$(diff \
		<(grep "distance =" correlation_global_chainx_$m     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
		<(grep "anchored edit distance" correlation_global_llchain_$m | cut -d' ' -f16))
	if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
done
echo " done (no difference)."

# time, memory, iterations
echo "similarity" > stats_correlation_headers

for t in correlation_global_chainx correlation_global_chainx_opt correlation_global_llchain
do
	echo "time (s)" > stats_correlation_time_$t
	echo "space (kb)" > stats_correlation_space_$t
done

for m in "90_100" "80_90" "75_80"
do
	echo "$m" >> stats_correlation_headers
	for t in correlation_global_chainx correlation_global_chainx_opt correlation_global_llchain
	do
		grep "total time elapsed" ${t}_${m} | cut -d' ' -f1 >> stats_correlation_time_$t
		grep "maxresident" ${t}_${m} | cut -d' ' -f1 >> stats_correlation_space_$t
	done
done
paste -d'$' stats_correlation_headers \
	stats_correlation_time_correlation_global_chainx     stats_correlation_space_correlation_global_chainx \
	stats_correlation_time_correlation_global_chainx_opt stats_correlation_space_correlation_global_chainx_opt \
	stats_correlation_time_correlation_global_llchain stats_correlation_space_correlation_global_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

echo "# Semiglobal sequence comparison, 100 queries"
usrbintimeoptions="%e total time elapsed (s)\n%M maxresident k"
ref="$thisfolder/../../ext/ChainX/data/correlation_semiglobal/e_coli_DH1_reference.fasta"
llchainmode="--mode semiglobal"
chainxmode="-m sg"

# run ChainX*
for m in "75_80" "80_90" "90_100"
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/correlation_semiglobal/count100_mutated_$m.fa" \
		>> correlation_semiglobal_chainx_$m 2>> correlation_semiglobal_chainx_$m
done

# run ChainX-opt*
for m in "75_80" "80_90" "90_100"
do
	$usrbintime "$usrbintimeoptions" $chainx $chainxmode -t $ref -q "$thisfolder/../../ext/ChainX/data/correlation_semiglobal/count100_mutated_$m.fa" \
		--optimal >> correlation_semiglobal_chainx_opt_$m 2>> correlation_semiglobal_chainx_opt_$m
done

# run llchain
for m in "75_80" "80_90" "90_100"
do
	$usrbintime "$usrbintimeoptions" $llchain $llchainmode -t $ref -q "$thisfolder/../../ext/ChainX/data/correlation_semiglobal/count100_mutated_$m.fa" \
		>> correlation_semiglobal_llchain_$m 2>> correlation_semiglobal_llchain_$m
done

echo -n "Checking if the optimal chaining cost differs..."
for m in "75_80" "80_90" "90_100"
do
	check=$(diff \
		<(grep "distance =" correlation_semiglobal_chainx_$m     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
		<(grep "distance =" correlation_semiglobal_chainx_opt_$m | cut -d'=' -f2 | cut -d')' -f2 | tr -d " "))
	if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
	check=$(diff \
		<(grep "distance =" correlation_semiglobal_chainx_$m     | cut -d'=' -f2 | cut -d')' -f2 | tr -d " ") \
		<(grep "anchored edit distance" correlation_semiglobal_llchain_$m | cut -d' ' -f16))
	if [ "$check" != "" ] ; then echo " it differs!" ; exit 1 ; fi
done
echo " done (no difference)."

# time, memory, iterations
echo "similarity" > stats_correlation_headers

for t in correlation_semiglobal_chainx correlation_semiglobal_chainx_opt correlation_semiglobal_llchain
do
	echo "time (s)" > stats_correlation_time_$t
	echo "space (kb)" > stats_correlation_space_$t
done

for m in "90_100" "80_90" "75_80"
do
	echo "$m" >> stats_correlation_headers
	for t in correlation_semiglobal_chainx correlation_semiglobal_chainx_opt correlation_semiglobal_llchain
	do
		grep "total time elapsed" ${t}_${m} | cut -d' ' -f1 >> stats_correlation_time_$t
		grep "maxresident" ${t}_${m} | cut -d' ' -f1 >> stats_correlation_space_$t
	done
done
paste -d'$' stats_correlation_headers \
	stats_correlation_time_correlation_semiglobal_chainx     stats_correlation_space_correlation_semiglobal_chainx \
	stats_correlation_time_correlation_semiglobal_chainx_opt stats_correlation_space_correlation_semiglobal_chainx_opt \
	stats_correlation_time_correlation_semiglobal_llchain stats_correlation_space_correlation_semiglobal_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'
