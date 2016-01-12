#!/bin/bash


# go to directory
pushd ../../tutorials/2-creating_microkinetic_model 

# run test and grep out critical string
#for fn in CO_oxidation.mkm CO_oxidation2.mkm CO_oxidation_square.mkm
for fn in CO_oxidation2.mkm
do
    seed=$(basename ${fn} .mkm)
    outdir="translated_${seed}_local_smart"

    catmap to_kmc ${fn}
    kmos export translated_${seed}.xml -o

    cp ${fn} ${outdir}
    pushd ${outdir}

    #catmap run_kmc
    pwd
    ls


    popd

done

popd
