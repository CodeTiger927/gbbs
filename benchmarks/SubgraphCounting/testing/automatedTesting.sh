#!/bin/bash
#After putting in the info, use the command
#./automatedTesting.sh|sed '/^#/d'|sed '/^D/d'|sed '/^T/d'|sed '/^$/d'>./results/res_GRAPHNAME.txt
#Change GRAPHNAME to the graph name(like DBLP)
#If you don't have the results folder, please use the command "mkdir results"

fname="../inputs/graph_test_DBLP.txt"
types=(0)
nodes=500000
parameter2=(3)
batchesSize=(10 20 50 100 200)
useP=(0)

for a in ${types[@]}; do
	for b in ${batchesSize[@]}; do
		for c in ${parameter2[@]}; do
			for d in ${useP[@]}; do
				../SubgraphCounting -s -rounds 1 -type $a -size $b -edges $c -useP $d -scriptMode 1 $fname
			done
		done
	done
done
