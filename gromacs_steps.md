# MD Simulation for protein and ligand intraction using gromacs

Step 1: Protein and ligand with chimerax or with other software.

Step 2: Protein topology and ligand topology preparation
* 1. Prepare the protein topology with pdb2gmx
* 2. Prepare the ligand topology using external tools


#### Prepare ligand topology with swissparam or Cgenff. We get 7 output file from Swissparam.

 > Open LIG.mol2 file on the header bellow line After the ```@<TRIPOS>MOLECULE``` there is some name like ***.pdb or **** or any other. Replace that name with LIG and now file should look like.
`@<TRIPOS>MOLECULE`\
LIG



This commond fix the bond order in .mol2 file which is neccessory to generate the topoplogy. The output file LIG.mol2.
> 'perl sort_mol2_bonds.pl LIG.mol2 LIG.mol2'

> download the sort_mol2_bonds.pl from gromacs\
keep the .pl file in same folder with .mol2 file then use 'perl sort_mol2_bonds.pl LIG.mol2 LIG.mol2'

Now upload LIG.mol2 file to Swissparam It would generate 7 output file.

convert .pdb to gromacs format file .gro

`gmx editconf -f LIG.pdb -o LIG.gro`


#### Protein topology preparation

`gmx pdb2gmx -f REC.pdb -ignh`

Press the corresponding number to the force field and water model you like to choose 
this will generate toptopol.top, conf.gro, posre.itp

Step 3:
#### Build the Complex
> The conf.gro file contains the processed, force field-compliant structure of our protein, and LIG.gro contains the processed structure of the ligand. To build complex we need to copy the LIG.gro file content from 3rd line to second last line and paste it to the second last line of the conf.gro file.\
In conf.gro file there is a number which we need to update:
Number: Total line in after edit - 3
So cross check if complex built correctly we can open the file in chimerax. if it opens and we see the complex that means file is good.



gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral



gmx make_ndx -f lig.gro -o index_lig.ndx

gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000


gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr

gmx mdrun -deffnm nvt

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr

gmx mdrun -deffnm npt


gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0.1.tpr