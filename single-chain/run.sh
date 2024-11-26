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
        sub ${ion} ${P} ${A} ${U} ${mol}
        sleep .5 
    done
    }

#ions=(100 200 400)
ions=(20 100 200 400 600)
T=293

name=polyR30


main -80 136 ${name} &

#main -10 120 ${name} & 
#main -10 122 ${name} & 
#main -10 124 ${name} &
#main -10 126 ${name} &
#main -10 128 ${name} &

#main 0 112 ${name} &
#main 0 114 ${name} &
#main 0 116 ${name} &
#main 0 118 ${name} & 
#main 0 120 ${name} & 
#main 0 122 ${name} & 
#main 0 124 ${name} & 
#main 0 126 ${name} & 
#main 0 128 ${name} & 

#main 10 112 ${name} & 
#main 10 114 ${name} & 
#main 10 116 ${name} & 
#main 10 118 ${name} & 
#main 10 120 ${name} & 
#main 10 122 ${name} & 
#
#main 20 112 ${name} &
#main 20 114 ${name} &
#main 20 116 ${name} &
#main 20 118 ${name} &
#main 20 120 ${name} &
#main 20 122 ${name} &

#main 30 118 ${name} &
#main 30 120 ${name} &
#main 30 122 ${name} &

#main 40 108 ${name} &
#main 40 110 ${name} &
#main 40 112 ${name} &

#main 50 108 ${name} &
#main 50 110 ${name} &
#main 50 112 ${name} &

#main 60 108 ${name} &
#main 60 110 ${name} &
#main 60 112 ${name} &
