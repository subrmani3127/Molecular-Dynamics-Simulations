gmx2 pdb2gmx -f KPC2-78.pdb -o KPC2-78_processed.gro -water tip3p
gmx2 editconf -f KPC2-78_processed.gro -o KPC2-78_newbox.gro -c -d 1.5 -bt dodecahedron
gmx2 solvate -cp KPC2-78_newbox.gro -cs spc216.gro -o KPC2-78_solv.gro -p topol.top
gmx2 grompp -f ions.mdp -c KPC2-78_solv.gro -p topol.top -o ions_78.tpr
gmx2 genion -s ions_78.tpr -o KPC2-78_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx2 grompp -f minim.mdp -c KPC2-78_ions.gro -p topol.top -o em_78.tpr
gmx2 mdrun -v -deffnm em_78
gmx2 grompp -f nvt.mdp -c em_78.gro -r em_78.gro -p topol.top -o KPC2-78_nvt.tpr
gmx2 mdrun -v -deffnm KPC2-78_nvt
gmx2 grompp -f npt.mdp -c KPC2-78_nvt.gro -r KPC2-78_nvt.gro -t KPC2-78_nvt.cpt -p topol.top -o KPC2-78_npt.tpr
gmx2 mdrun -v -deffnm KPC2-78_npt
gmx2 grompp -f md.mdp -c KPC2-78_npt.gro -t KPC2-78_npt.cpt -p topol.top -o KPC2-78_MD.tpr
gmx2 mdrun -v -deffnm KPC2-78_MD
