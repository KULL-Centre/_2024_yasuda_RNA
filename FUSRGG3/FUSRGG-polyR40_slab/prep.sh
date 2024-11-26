#! /bin/bash


function main(){
	protein=$1
	rna=$2
	box=$3
    mkdir  FUSRGG-${protein}_polyR40-${rna}
	cd FUSRGG-${protein}_polyR40-${rna}
	bash ../src/write_prepare_full.sh  ${protein} ${rna} ${box}
    cp ../src/fastabib.fasta ./fastabib.fasta 
    cp ../src/residues.csv ./residues.csv 
	python3 prepare_full.py
	cd ../
    #rm -r ${mol}_${hp}_${ion} #for test
	}

#no llps
#main 500 5 [20,20,150]  
#main 500 10 [20,20,150]  
#main 500 20 [20,20,150]  

#llps
main 500 40 [20,20,150]  
main 500 80 [20,20,150]  
main 500 100 [20,20,150]  
main 500 120 [20,20,150]  
main 500 130 [20,20,150]  
main 500 140 [20,20,150]  
main 500 160 [20,20,150]  
main 500 200 [20,20,150]  

#no llps
#main 500 210 [20,20,150]
#main 500 240 [20,20,150]  
#main 500 320 [20,20,150]  
#main 500 500 [20,20,150]  
#main 500 550 [20,20,150]  
