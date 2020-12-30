#!/bin/bash
#After putting in the info, use the command
#./automatedTesting.sh
#Change GRAPHNAME to the graph name(like DBLP)
#If you don't have the results folder, please use the command "mkdir results"

timestamp=`date +%y%m%d%H%M%S`
fname="../inputs/graph_test_DBLP.txt"
types=(0)
nodes=20000
parameter2=(5)
batchesSize=(10, 20, 30, 40, 50)
useP=(0, 1)

for a in ${types[@]}; do
	for b in ${batchesSize[@]}; do
		for c in ${parameter2[@]}; do
			for d in ${useP[@]}; do
<<<<<<< HEAD
				../SubgraphCounting -s -rounds 3 -type $a -size $b -edges $c -nodes $nodes -P $d -scriptMode 1 $fname |sed '/^#/d'|sed '/^D/d'|sed '/^T/d'|sed '/^$/d'>>./results/res_DBLP_$timestamp.txt
=======
				../SubgraphCounting -s -rounds 3 -type $a -size $b -edges $c -nodes $nodes -useP $d -scriptMode 1 $fname |sed '/^#/d'|sed '/^D/d'|sed '/^T/d'|sed '/^$/d'>>./results/res_DBLP.txt
>>>>>>> c06eee9132d58e990ea5e9a511b390c699bf4dc8
			done
		done
	done
done
