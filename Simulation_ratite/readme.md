The folder contains simulation data and result for bird phylogeny. **simu_500_200_diffr_\*.bed** and **simu_500_200_diffr_\*.fasta** are input data. 1-1 corresponds to scenario (1) (all branches conserved) in the main paper; 2-* corresponds to scenarios (2)~(9) (different acceleration patterns). In the bed file, the 1th column is element name (ID); the 2nd and 3rd columns are start and end positions in the alignment for each element; the 6th and 7th columns are conserved and accelerated rates to generate the DNA sequences for that element. [param2-6.txt](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation_ratite/param2-6.txt) is the parameter file to get the result for scenario (7) (all ratites' branches are accelerated). You could try to run: 
```bash
PhyloAcc Simulation_ratite/param2-6.txt
```
under the root directory and it will generate results in [result_phyloAcc](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation_ratite/result_phyloAcc). 

Results from PhyloAcc are shown in [result_phyloAcc/](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation_ratite/result_phyloAcc). Results from phyloP are in  [result_phyloP/](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation_ratite/result_phyloP). Input files and Results using PAML and phytools are in [PAML/](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation_ratite/PAML). Input files and Results using BEAST2 are in [BEAST/](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation_ratite/BEAST).

