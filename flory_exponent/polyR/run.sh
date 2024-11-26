#! /bin/bash


function sub(){
    ion=$1
    P=$2
    A=$3
    U=$4
    mol=$5

    hp=p${P}a${A}u${U}
    cd ${mol}_${hp}_${ion}/${mol}/${T}
    python3 run.py
    cd ../../../
    }



function main(){
    P=$1
    A=$2
    U=$2
    mol=$3

    for ion in ${ions[@]};do
        sub ${ion} ${P} ${A} ${U} ${mol} & 
        sleep .5 
    done
    }

#ions=(100 200 400)
ions=(225 525)
T=293

name=polyR10
main 0 118 ${name} &
name=polyR20
main 0 118 ${name} &
name=polyR30
main 0 118 ${name} & 
name=polyR40
main 0 118 ${name} &
name=polyR50
main 0 118 ${name} &
name=polyR60
main 0 118 ${name} &
name=polyR70
main 0 118 ${name} &
name=polyR80
main 0 118 ${name} &
name=polyR90
main 0 118 ${name} &
name=polyR100
main 0 118 ${name} &
