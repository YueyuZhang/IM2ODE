#!/bin/bash
#submit all the PWmat jobs and queue up

i=0
n=4 #number of GPU

for x in `ls POSCAR*`
do 
echo $x


	ind=`echo $x | sed 's/POSCAR//'`
        echo run$ind
	rm -rf  pwmat_$ind
	mkdir pwmat_$ind
	cp etot.input *.UPF check.x run convert_to_config.x MV2CONTCAR.x writekp.py pwmat_$ind
        mv atom.config$ind pwmat_$ind/atom.config
	cp POSCAR$ind pwmat_$ind/POSCAR.vasp
	cp POSCAR$ind pwmat_$ind/POSCAR

	stat=0
        while [ $stat -eq 0 ]
        do
            sleep 2s
            stat=1
            ps -l | grep PWmat >back
            m=$(wc -l back| awk '{printf("%d\n", $1)}')
            echo m=$m
            rm back
            if [[ $m -gt $(($n-1)) ]]; then
                stat=0
            fi
        done

	cd pwmat_$ind
#        ./convert_to_config.x< POSCAR.vasp
        python writekp.py 0.08
        cat KPOINTS >> etot.input
        ./check.x >>tmp
	./run&

	let "i=$i+1"
        cd ..

done


	stat=0
        while [ $stat -eq 0 ]; do
            sleep 5s
            stat=1
            ps -l | grep PWmat >back
            m=$(wc -l back| awk '{printf("%d\n", $1)}')
            echo *m=$m
            rm back
            if [ $m -gt 0 ]; then
                stat=0
            fi
        done


