unzip charmm36-feb2021.ff
gmx pdb2gmx -f O08324.pdb -o O08324_processed.gro -water spce
gmx editconf -f O08324_processed.gro -o O08324_box.gro -c -d 1.5 -bt dodecahedron
gmx insert-molecules -f O08324_box.gro -ci bgc.gro -nmol 41 -try 20 -o O08324_bgc_box.gro
gmx2 solvate -cp O08324_bgcbox.gro -cs spc216.gro -o O08324_solv.gro -p topol.top
#vi topol.top  # at top bgc.pr m # bottom bgc.itp, number of  bgc
gmx2 grompp -f ions.mdp -c O08324_solv.gro -p topol.top -o ions.tpr

gmx2 genion -s ions.tpr -o O08324_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx2 grompp -f em.mdp -c O08324_ions.gro -p topol.top -o em.tpr
gmx2 mdrun -gpu_id 0 -v -deffnm em

gmx2 make_ndx -f em.gro -o index.ndx <<EOF
1 | 13
q
EOF
gmx2 energy -f em.edr <<EOF
11 0
EOF

gmx2 make_ndx -f bgc.gro -o index_bgc.ndx <<EOF
0 & ! a H*
q
EOF

gmx2 genrestr -f bgc.gro -n index_bgc.ndx -o posre_bgc.itp -fc 1000 1000 1000 <<EOF
1
EOF

gmx2 grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -n index.ndx
gmx2 mdrun -gpu_id 0 -v -deffnm nvt

gmx2 grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -n index.ndx
gmx2 mdrun -gpu_id 0 -v -deffnm npt
