#!/bin/bash
#submit all the vasp jobs
for x in `ls POSCAR*`
do
    ind=`echo $x | sed 's/POSCAR//'`
    rm -rf vasp_$ind
    mkdir vasp_$ind
    cp INCAR_* KPOINTS POTCAR vasp.pbs writekp.py run vasp_$ind
    cp POSCAR$ind vasp_$ind/POSCAR
    cp POSCAR$ind vasp_$ind/POSCAR.bak
    #For magnetic system
    if [ -s MAGMOM_$ind ]; then
       cp MAGMOM_$ind vasp_$ind/MAGMOM
       (
	   cd vasp_$ind
	   for y in `ls INCAR_*`
	   do
	       cat MAGMOM >> $y
	   done
       )
    fi
    cd vasp_$ind
 #   sbatch vasp.pbs |tail -1|awk '{print $NF}' > job_ID.dat
#PBS
    oldpwd=`pwd`
    ssh -p 2323 node1 "cd $oldpwd ; ./run |grep 'node2'|head -1 " > job_ID.dat
    sleep 2s
    cd ..
done

#now check whether all the calculations finish

stat=0
while [ $stat -eq 0 ]; do
    sleep 60s
    stat=1
    /opt/gridview/pbs/dispatcher/bin/qstat | grep xggong | awk '{printf("%6d\n", $1)}' >back
    n=$(wc -l back | awk '{printf("%d\n", $1)}')
    nu=$(($n+1))
    for x in `ls -d vasp_*/`
    do
        jobid=$(head -1 $x/job_ID.dat|awk '{printf("%6d\n", $1)}')
        for((l=1;l<$nu;l++)); do
            jobone=$(sed -n "${l}p" back | awk '{printf("%6d\n", $1)}')
            if [[ $jobid == $jobone ]]; then
                stat=0
            fi
        done
    done
    rm back
done

rm -rf POSCAR_*
