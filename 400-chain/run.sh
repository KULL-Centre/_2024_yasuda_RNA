#! /bin/bash

function sub(){
    ion=$1
    P=$2
    A=$3
    U=$4
    mol=$5

	hp=p${P}a${A}u${U}
    cd ${mol}_${hp}_${ion}/${mol}/${T}
    #python3 run.py 
    python3 post.py ${mol} ${n_chains}
    cd ../../../

	}		

function main(){
	P=$1
	A=$2
	U=$3
	mol=$4
    ion=$5
    sub ${ion} ${P} ${A} ${U} ${mol}
	}

n_chains=400
T=293
name=polyR30

#l1=0
#l2=118
#ion=20
#main ${l1} ${l2} ${l2} ${name} ${ion}
#ion=100
#main ${l1} ${l2} ${l2} ${name} ${ion}
#ion=200
#main ${l1} ${l2} ${l2} ${name} ${ion}
#ion=600
#main ${l1} ${l2} ${l2} ${name} ${ion}
#ion=400
#main ${l1} ${l2} ${l2} ${name} ${ion}


