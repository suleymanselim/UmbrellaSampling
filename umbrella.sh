#!/bin/bash
export gmx=/home/sulejmani/Downloads/Programs/gromacs-19_gpu/bin/gmx
echo 0 | $gmx trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep
for ps in `seq 1 500`;do 
$gmx distance -s pull.tpr -f conf${ps}.gro -n index.ndx -select 'com of group "Chain_L" plus com of group "Chain_R"' -oall dist${ps}.xvg 
d=`tail -n 1 dist${ps}.xvg | awk '{print $2}'`
echo "${ps} ${d}" >> summary_distances.dat
rm dist${ps}.xvg
done
range=`cat summary_distances.dat| awk '{if(min==""){min=max=$2}; if($2>max) {max=$2}; if($2<min) {min=$2}; total+=$2; count+=1} END {print min" 0.2 "max}'`
for xx in `seq -f "%f" ${range}`;do 
aa=`awk -v c=2 -v t=${xx} '{a[NR]=$c}END{
        asort(a);d=a[NR]-t;d=d<0?-d:d;v = a[NR]
        for(i=NR-1;i>=1;i--){
                m=a[i]-t;m=m<0?-m:m
                if(m<d){
                    d=m;v=a[i]
                }
        }
        print v
}' summary_distances.dat`
grep $aa summary_distances.dat |head -n1 >> list
done
awk '!seen[$0]++' list
aa=1
while read i;do num=`echo $i|awk '{print $1}'`;mv conf${num}.gro umb${aa}.gro;aa=`echo $aa|awk '{print $1 + 1}'`;done <list
rm conf*

for ii in `seq 1 26`;do
$gmx grompp -f files/npt_umbrella.mdp -c files/umb${ii}.gro -p topol.top -r conf${ii}.gro -n index.ndx -o npt${ii}.tpr
gmx mdrun -deffnm npt${ii}
gmx grompp -f files/md_umbrella.mdp -c npt${ii}.gro -t npt${ii}.cpt -p topol.top -r npt${ii}.gro -n index.ndx -o umbrella${ii}.tpr
gmx mdrun -deffnm umbrella${ii}
done
