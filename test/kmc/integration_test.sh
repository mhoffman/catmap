#!/bin/bash -eux


# go to directory
pushd ../../tutorials/kmc_sandbox

# run test and grep out critical string
#for fn in CO_oxidation.mkm CO_oxidation2.mkm CO_oxidation_square.mkm

#for fn in CO_oxidation.mkm CO_oxidation2.mkm CO_oxidation3.mkm
for fn in CO_oxidation.mkm CO_oxidation2.mkm
do
    seed=$(basename ${fn} .mkm)
    outdir="${seed}_kmc_local_smart"

    catmap to_kmc ${fn}
    
    if ! diff -q ${seed}_kmc.ini reference/${seed}_kmc.ini > /dev/null
    then
        echo -e "\n\nGenerated kmc file for ${seed} (${fn}) differs from reference."
        diff ${seed}_kmc.ini reference/${seed}_kmc.ini
        exit 1
    fi


    kmos export -s ${seed}_kmc.ini -o

    cp ${fn} ${outdir}
    pushd ${outdir}

    #catmap run_kmc
    pwd
    ls -rtlh


    popd

done

popd

echo "kMC integration test finished successful."
