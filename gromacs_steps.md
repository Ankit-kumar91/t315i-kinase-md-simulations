# MD Simulation for protein and ligand intraction using gromacs

**Step 1:** Protein and ligand with chimerax or with other software.

**Step 2:** Protein topology and ligand topology preparation
* 1. Prepare the protein topology with pdb2gmx
* 2. Prepare the ligand topology using external tools


#### Prepare ligand topology with swissparam or Cgenff. We get 7 output file from Swissparam.

 > Open LIG.mol2 file on the header bellow line After the ```@<TRIPOS>MOLECULE``` there is some name like ***.pdb or **** or any other. Replace that name with LIG and now file should look like.

```
@<TRIPOS>MOLECULE
LIG
```

This commond fix the bond order in .mol2 file which is neccessory to generate the topoplogy. The output file LIG.mol2.
> 'perl sort_mol2_bonds.pl LIG.mol2 LIG.mol2'

> download the sort_mol2_bonds.pl from gromacs\
keep the .pl file in same folder with .mol2 file then use 'perl sort_mol2_bonds.pl LIG.mol2 LIG.mol2'

Now upload LIG.mol2 file to Swissparam It would generate 7 output file.

convert .pdb to gromacs format file .gro

`gmx editconf -f lig.pdb -o lig.gro`


#### Protein topology preparation

`gmx pdb2gmx -f REC.pdb -ignh`

Press the corresponding number to the force field and water model you like to choose 
this will generate toptopol.top, conf.gro, posre.itp

**Step 3:**
#### Build the Complex
> The conf.gro/complex.gro file contains the processed, force field-compliant structure of our protein, and LIG.gro contains the processed structure of the ligand. To build complex we need to copy the LIG.gro file content from 3rd line to second last line and paste it to the second last line of the conf.gro file.\
In conf.gro file there is a number which needs to be update:\
Number: Total line in after edit - 3\
To confirm that, the complex built correctly we can open the file in chimerax. if it opens and we see the complex that means complex built correctly.

* Now we need to do some changes in topol.top file.
> bellow 
```
; Include forcefield parameters
#include "charmm27.ff/forcefield.itp"
```

> we need to add\
```
; Include ligand parameters
#include "lig.itp
```

> At the bottom of the file add \
LIG        1 

Edit the lig.itp

bellow these line

```
[ moleculetype ]
; Name nrexcl 
add LIG 3
```

if LIG 3 is already present then no need to edit.

**Step 4:**
## Box creation

> gmx editconf -f conf.gro -o box.gro -bt dodecahedron -d 1.2

## Solvation

> gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o sol.gro

**Step 5:**
## Add ions

> gmx grompp -f ions.mdp -c sol.gro -p topol.top -o ions.tpr

if you want maintain perticular conc use this command. 
> gmx genion -s ions.tpr -o sol_ions.gro -p topol.top -pname NA -nname CL -conc 0.15 -neutral

else
> gmx genion -s ions.tpr -o sol_ions.gro -p topol.top -pname NA -nname CL -neutral

Select a group: 15

Enter number correspond to SOL 


**Step 6:**
## Energy Minimization

> gmx grompp -f em.mdp -c sol_ions.gro -p topol.top -o em.tpr

> gmx mdrun -v -deffnm em

**Step 6:**
## Equilibration


Make index file for ligand
> gmx make_ndx -f lig.gro -o index_lig.ndx

then type
> 0 & ! a H*

then press
> q

 - Restraining the Ligand

> gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000

select group number 3
> (System_&_!H)

* Open topol.top 

After 
```cpp
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif
```

add this 
```
; Ligand position restraints
#ifdef POSRES
#include "posre_lig.itp"
#endif
```
* Make index file for system

> gmx make_ndx -f em.gro -o index.ndx

then type
> 1 | 13 (Protein | LIG)

then press
> q

**Step 6:**
## NVT equilibration

> gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr

> gmx mdrun -deffnm nvt

gpu
gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu

## NPT equilibration
> gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr

> gmx mdrun -deffnm npt

for gpu 
gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu

**Step 6:** 
## Production run
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0.1.tpr

mx mdrun -deffnm md -nb gpu -pme gpu -bonded gpu