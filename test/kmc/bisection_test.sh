#!/bin/bash


# go to directory
pushd ../../tutorials/2-creating_microkinetic_model 

# run test and grep out critical string
for fn in CO_oxidation.mkm CO_oxidation2.mkm CO_oxidation_square.mkm
do
    catmap to_kmc ${fn} | grep DLASCL
    if [ "$?" -eq "0" ]; then
         echo "Text found."
        popd 
        exit 1
    else
         echo "Text not found."
        popd 
        exit 0
    fi
done

