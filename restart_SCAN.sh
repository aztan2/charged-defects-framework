cd $1
for q in '1'; do
    mkdir -p charge_$q
    cd charge_$q
    for cell in '4x4x1'; do
        mkdir -p $cell
        cd $cell
        for vac in '15' '20'; do
            mkdir -p vac_$vac
            cd vac_$vac
            echo $q $cell $vac
            cp ../../../../../GGA/mag/charge_$q/$cell/vac_$vac/{INCAR,KPOINTS,POTCAR,CONTCAR,WAVECAR,sub*,def*} .
            mv CONTCAR POSCAR
            python /ufrc/hennig/annemarietan/scripts/gen_incar.py $q --runtype relax --functional SCAN+rVV10
            python /ufrc/hennig/annemarietan/scripts/gen_submit.py --nodes 2 --time 4-00:00:00 
            if [ -d ../../../../../GGA/mag/charge_$q/$cell/vac_$vac/bulkref ]; then
                mkdir -p bulkref
                cd bulkref
                cp ../../../../../../GGA/mag/charge_$q/$cell/vac_$vac/bulkref/{INCAR,KPOINTS,POTCAR,CONTCAR,WAVECAR,sub*} .
                mv CONTCAR POSCAR
                python /ufrc/hennig/annemarietan/scripts/gen_incar.py $q --runtype relax --functional SCAN+rVV10
                cd ../
            fi
            cd ../
        done
        cd ../
    done
    cd ../
done
