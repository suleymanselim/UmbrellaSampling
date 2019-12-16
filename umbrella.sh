#!/bin/bash
export gmx=~/Programs/gromacs/bin/gmx           #Gromacs directory
export files=~/files                            #Directory containing input files

#Generating a molecular topology for an umbrella sampling simulation
echo 1|$gmx pdb2gmx -f $files/complex.pdb -o complex.gro -water tip3p 

$gmx editconf -f complex.gro -o newbox.gro -center 2.5 3 3 -box 15 6 6
$gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
$gmx grompp -f $files/ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL|$gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
$gmx grompp -f $files/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
$gmx mdrun -deffnm em

gmx grompp -f $files/npt.mdp -c em.gro -p topol.top -r em.gro -o npt.tpr
gmx mdrun -deffnm npt

$gmx grompp -f $files/md_pull.mdp -c npt.gro -p topol.top -r npt.gro -n index.ndx -t npt.cpt -o pull.tpr
$gmx mdrun -deffnm pull -pf pullf.xvg -px pullx.xvg

echo 0 | $gmx trjconv -s pull.tpr -f pull.xtc -o conf.gro -sep

for ps in `seq 1 500`;do 
$gmx distance -s pull.tpr -f conf${ps}.gro -n index.ndx -select 'com of group "Chain_L" plus com of group "Chain_R"' -oall dist${ps}.xvg 
d=`tail -n 1 dist${ps}.xvg | awk '{print $2}'`
echo "${ps} ${d}" >> distances.dat
rm dist${ps}.xvg
done

range=`cat distances.dat| awk '{if(min==""){min=max=$2}; if($2>max) {max=$2}; if($2<min) {min=$2}; total+=$2; count+=1} END {print min" 0.2 "max}'`

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
}' distances.dat`
grep $aa distances.dat |head -n1 >> list
done

aa=1
awk '!a[$0]++' list|while read i;do num=`echo $i|awk '{print $1}'`;mv conf${num}.gro $files/umb${aa}.gro;aa=`expr $aa+ 1`;done
rm conf*

for ii in $(seq 1 `awk '!a[$0]++' list |wc -l`);do
$gmx grompp -f $files/npt_umbrella.mdp -c $files/umb${ii}.gro -p topol.top -r $files/umb${ii}.gro -n index.ndx -o npt${ii}.tpr
$gmx mdrun -deffnm npt${ii}
$gmx grompp -f $files/md_umbrella.mdp -c npt${ii}.gro -t npt${ii}.cpt -p topol.top -r npt${ii}.gro -n index.ndx -o umbrella${ii}.tpr
$gmx mdrun -deffnm umbrella${ii}
done
