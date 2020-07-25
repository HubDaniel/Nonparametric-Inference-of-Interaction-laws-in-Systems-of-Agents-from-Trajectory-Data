## Nonparametric inference of interaction laws in systems of agents from trajectory data
## Fei Lu, Ming Zhong, Sui Tang, and Mauro Maggioni

Reproducing the result in R. 
Original Code in Matlab: https://github.com/MingZhongCodes    (open source)

All code in R is reproduced by Xiaoran Ma with necessary adjustments.

Any suggestions will be appreciated.

This the 2nd Version:
1. Version 2 can be used for heterogeneous cases
2. Reformat Version 2 such that homogeneous cases can also be run
3. Cross check (without obs. noise) the result for homogeneous case use Version 1 and Version 2 (scale quantity function error fixed)
4. Reformat Version 2 to use packaged initial information
5. Add trajectory plots and predicted trajectory plots
6. Fixed trajectory error calculation problem
7. Store function in nested list for heterogeneous cases. Deal with lazy evaluation problem (lapply for simple list and force function for nested list)
8. Test homogeneous case with sine function
9. Add tryCatch to ode function in case cannot be solved
10. Use ode45 will be faster in generating trajectories, in general
11. Add estimated rhoLT plot
