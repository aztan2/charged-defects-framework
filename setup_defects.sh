cd $1
for q in '0' '1' '-1' '-2'; do
    mkdir -p charge_$q
    cd charge_$q
    for cell in '3x3x1' '4x4x1'; do
        mkdir -p $cell
        cd $cell
        for vac in '15' '20'; do
            mkdir -p vac_$vac
            cd vac_$vac
            echo $q $cell $vac
            python /ufrc/hennig/annemarietan/scripts/gen_defect_supercells.py $1 $q $cell $vac
            cp $1/POTCAR .
            python /ufrc/hennig/annemarietan/scripts/gen_incar.py $q --runtype relax --functional PBE
            python /ufrc/hennig/annemarietan/scripts/gen_kpts_grid.py --kppa 440
            python /ufrc/hennig/annemarietan/scripts/gen_submit.py --nodes 1 --time 24:00:00
            if [ -d $1/charge_$q/$cell/vac_$vac/bulkref ]; then
                cd bulkref
                cp $1/POTCAR .
                python /ufrc/hennig/annemarietan/scripts/gen_incar.py $q
                cp ../{KPOINTS,submitVASP.sh} .
                cd ../
            fi 
        cd ../
        done
        cd ../
    done
    cd ../
done
