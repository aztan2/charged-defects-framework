cd $1
for q in '-1' '1'; do
    cd charge_$q
    for cell in '3x3x1' '4x4x1'; do
        cd $cell
        for vac in '15' '20'; do
            echo $q $cell $vac
            cd vac_$vac
#            cd vac_$vac/dos
            python /ufrc/hennig/annemarietan/scripts/gen_SPHInX_input_file.py 17.18 5.40
#            python /ufrc/hennig/annemarietan/scripts/gen_SPHInX_input_file.py 16.74 5.93    
#            python /ufrc/hennig/annemarietan/scripts/gen_SPHInX_input_file.py 16.04 5.88
            dir_locpot=$(pwd)
            cd correction
            python /ufrc/hennig/annemarietan/scripts/get_alignment_correction_2d.py $2/charge_0/$cell/vac_$vac/bulkref/LOCPOT $dir_locpot/LOCPOT 520 $q --logfile getalign.log
#            python /ufrc/hennig/annemarietan/scripts/get_alignment_correction_2d.py $2/charge_0/$cell/vac_$vac/bulkref/dos/LOCPOT $dir_locpot/LOCPOT 520 $q --logfile getalign.log
            cd ../../
#            cd ../../../
        done        
        cd ../
    done
    cd ../
done
