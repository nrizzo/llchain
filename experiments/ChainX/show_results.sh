#!/bin/bash
set -euo pipefail
thisfolder=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # https://stackoverflow.com/questions/59895/how-do-i-get-the-directory-where-a-bash-script-is-located-from-within-the-script
cd $thisfolder/output

echo "# Global sequence comparison"
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
paste -d'$' stats_time_global_mums \
	stats_time_time_global_chainx             stats_space_time_global_chainx \
	stats_time_time_global_chainx_opt         stats_space_time_global_chainx_opt  \
	stats_time_time_global_llchain            stats_space_time_global_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'
echo "# Semiglobal sequence comparison"
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
paste -d'$' stats_time_semiglobal_mums \
	stats_time_time_semiglobal_chainx             stats_space_time_semiglobal_chainx \
	stats_time_time_semiglobal_chainx_opt         stats_space_time_semiglobal_chainx_opt  \
	stats_time_time_semiglobal_llchain            stats_space_time_semiglobal_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'
echo "# Global sequence comparison, 100 queries"
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
paste -d'$' stats_correlation_headers \
	stats_correlation_time_correlation_global_chainx     stats_correlation_space_correlation_global_chainx \
	stats_correlation_time_correlation_global_chainx_opt stats_correlation_space_correlation_global_chainx_opt \
	stats_correlation_time_correlation_global_llchain stats_correlation_space_correlation_global_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'

echo "# Semiglobal sequence comparison, 100 queries"
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
paste -d'$' stats_correlation_headers \
	stats_correlation_time_correlation_semiglobal_chainx     stats_correlation_space_correlation_semiglobal_chainx \
	stats_correlation_time_correlation_semiglobal_chainx_opt stats_correlation_space_correlation_semiglobal_chainx_opt \
	stats_correlation_time_correlation_semiglobal_llchain stats_correlation_space_correlation_semiglobal_llchain \
	| cat <(echo -e "\$ChainX*\$\$ChainX-opt*\$\$llchain") - | column -t -s'$'
