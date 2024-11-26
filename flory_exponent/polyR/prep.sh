#! /bin/bash

# simulation conditions
ions=(225 525) 
#Ps=(0 1 2 3 4 5 6 7 8 9 10)


function main(){
	ion=$1
	P=$2
	A=$3
	U=$3
	mol=$4

	hp=p${P}a${A}u${U}

	mkdir -p ${mol}_${hp}_${ion}

	cd ${mol}_${hp}_${ion}

	# make residue file
	p=`echo "scale=5; ${P} / 100 " | bc`
	a=`echo "scale=5; ${A} / 100 " | bc`
	u=`echo "scale=5; ${U} / 100 " | bc`
	bash ../src/write_res.sh ${p} ${a} ${u}

	# make prepare file
	bash ../src/write_prepare_full.sh ${mol} ${hp} ${ion}
    
    # preparer fasta
    cp ../src/fasta.fasta ./fastabib.fasta 
 
    # run prepare file
	python3 prepare_full.py
	cd ../
    #rm -r ${mol}_${hp}_${ion} #for test
    echo ${hp} ${ion}  
	}

for ion in ${ions[@]};do
name=polyR10
main ${ion} 0 118 ${name}
name=polyR20
main ${ion} 0 118 ${name}
name=polyR30
main ${ion} 0 118 ${name}
name=polyR40
main ${ion} 0 118 ${name}
name=polyR50
main ${ion} 0 118 ${name}
name=polyR60
main ${ion} 0 118 ${name}
name=polyR70
main ${ion} 0 118 ${name}
name=polyR80
main ${ion} 0 118 ${name}
name=polyR90
main ${ion} 0 118 ${name}
name=polyR100
main ${ion} 0 118 ${name}
done

