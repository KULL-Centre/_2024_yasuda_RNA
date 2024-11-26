#! /bin/bash

# simulation conditions
ions=(20 100 200 400 600) 
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

name=polyR30
for ion in ${ions[@]};do
        main ${ion} -80 136 ${name}

#        main ${ion} -10 120 ${name} 
#        main ${ion} -10 122 ${name}
#        main ${ion} -10 124 ${name}
#        main ${ion} -10 126 ${name}
#        main ${ion} -10 128 ${name}


#        main ${ion} 0 112 ${name} 
#        main ${ion} 0 114 ${name}
#        main ${ion} 0 116 ${name}
#        main ${ion} 0 118 ${name} 
#        main ${ion} 0 120 ${name} 
#        main ${ion} 0 122 ${name}
#        main ${ion} 0 124 ${name}

#        main ${ion} 10 112 ${name} 
#        main ${ion} 10 114 ${name} 
#        main ${ion} 10 116 ${name} 
#        main ${ion} 10 118 ${name} 
#        main ${ion} 10 120 ${name} 
#        main ${ion} 10 122 ${name}  

#        main ${ion} 20 112 ${name} 
#        main ${ion} 20 114 ${name}
#        main ${ion} 20 116 ${name}
#        main ${ion} 20 118 ${name} 
#        main ${ion} 20 120 ${name} 
#        main ${ion} 20 122 ${name} 
        
#        main ${ion} 30 118 ${name} 
#        main ${ion} 30 120 ${name} 
#        main ${ion} 30 122 ${name} 
#        
#        main ${ion} 40 108 ${name} 
#        main ${ion} 40 110 ${name} 
#        main ${ion} 40 112 ${name}
#
#        main ${ion} 50 108 ${name}
#        main ${ion} 50 110 ${name}
#        main ${ion} 50 112 ${name}
#
#        main ${ion} 60 108 ${name}
#        main ${ion} 60 110 ${name}
#        main ${ion} 60 112 ${name}

 
done

