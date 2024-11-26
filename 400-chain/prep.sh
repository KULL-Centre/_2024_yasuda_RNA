#! /bin/bash


function main(){
	ion=$1
	P=$2
	A=$3
	U=$4
	mol=$5
    gpu=$6
    T=293

	hp=p${P}a${A}u${U}

	mkdir -p ${mol}_${hp}_${ion}

	cd ${mol}_${hp}_${ion}

	# make residue file
	p=`echo "scale=5; ${P} / 100 " | bc`
	a=`echo "scale=5; ${A} / 100 " | bc`
	u=`echo "scale=5; ${U} / 100 " | bc`
	bash ../src/write_res.sh ${p} ${a} ${u}

	# make prepare file
	bash ../src/write_prepare_full.sh ${mol} ${hp} ${ion} ${gpu}
    
    # preparer fasta
    cp ../src/fasta.fasta  ./fastabib.fasta

 
    # run prepare file
	python3 prepare_full.py
    
    # copy rdf calculation file (must be after runing prepare_full)
    cp ../src/analyse.py  ${mol}/${T}/
    cp ../src/post.py     ${mol}/${T}/
    p ../src/rdf_calc.py ${mol}/${T}/


	cd ../
    #rm -r ${mol}_${hp}_${ion} #for test
	}

name=polyR30

l1=0
l2=118    
ion=600
main ${ion} ${l1} ${l2} ${l2} ${name} 1
ion=400
main ${ion} ${l1} ${l2} ${l2} ${name} 1
ion=200
main ${ion} ${l1} ${l2} ${l2} ${name} 1
ion=100
main ${ion} ${l1} ${l2} ${l2} ${name} 1
ion=20
main ${ion} ${l1} ${l2} ${l2} ${name} 1
