The folder contains simulation data and example result. **simu_500_200_diffr_\*.bed** and **simu_500_200_diffr_\*.fasta** are input data. 1-1 corresponds to scenario (1) (all branches conserved) in the main paper; 2-* corresponds to scenarios (2)~(9) (different acceleration patterns). In the bed file, the 1th column is element name (ID); the 2nd and 3rd columns are start and end positions in the alignment for each element; the 7th and 8th columns are conserved and accelerated rates to generate the DNA sequences for that element. [param2-6.txt](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation/param2-6.txt) is the parameter file to get the result for scenario (7) (all ratites' branches are accelerated). You could try to run: 
```bash
PhyloAcc Simulation/param2-6.txt
```
under the root directory and example results are shown in [result](https://github.com/xyz111131/PhyloAcc/tree/master/Simulation/result). 
