cd $1
for q in '0' '1' '-1' '-2'; do
    cd charge_$q
    for cell in '3x3x1' '4x4x1'; do
        cd $cell
        for vac in '15' '20'; do
            cd vac_$vac
            echo $q $cell $vac
	    mkdir -p dos
	    cd dos
#            cp ../{INCAR,KPOINTS,POTCAR,CONTCAR,sub*,def*} .
            cp ../{INCAR,KPOINTS,POTCAR,CONTCAR,CHGCAR,sub*,def*} .
            mv CONTCAR POSCAR
            python /ufrc/hennig/annemarietan/scripts/gen_incar.py $q --runtype dos --functional PBE --soc
            python /ufrc/hennig/annemarietan/scripts/gen_submit.py --time 24:00:00 --vasp ncl_noz
	    cd ../
            if [ -d $1/charge_$q/$cell/vac_$vac/bulkref ]; then
                cd bulkref
	        mkdir -p dos
                cd dos
#                cp ../{INCAR,KPOINTS,POTCAR,CONTCAR,sub*} .
                cp ../{INCAR,KPOINTS,POTCAR,CONTCAR,CHGCAR,sub*,def*} .
                mv CONTCAR POSCAR
                python /ufrc/hennig/annemarietan/scripts/gen_incar.py $q --runtype dos --functional PBE --soc
                cd ../../
            fi
            cd ../
        done
        cd ../
    done
    cd ../
done
