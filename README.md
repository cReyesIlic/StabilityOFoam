# StabilityOFoam 1.0

## Solvers 
Solvers and cases to run a Jacobian-free matrix linear stability analysis.
Tested in OpenFOAM 7.x with cylinder wake (stable and unstable cases)

Steps for compiling
1. Install OpenFOAM (visit https://openfoam.org/)
2. Compiling icoFOAMSFD solver
3. Compiling icoFOAMLNS solver

### icoFOAMSFD solver

Solver icoFOAM based to calculate the base flow for stability analysis using selective frequency damping (SFD). It is necessary to define in transportProperties dictionary the valu of $\chi$ and rDelta (parameter of SFD).

### icoFOAMLNS

solver icoFOAM based which calculate the perturbation velocity by Linearized Navier Stokes,

## Linear Stability Analysis - Jacobian Free Method

To calculate the H matrix (size m x m, with m and user-specified parameter) it is necessary to follow the next steps:

1. Calculate the initial random $\zeta$ value, where $\zeta = [U_1\quad V_1 \quad W_1 \quad U_2 \quad V_2 \quad W_2 \quad ...]^T $
Begin Iteration k = 0:
  2. Calculate the new perturbation $U_{k + 1} = \beta\zeta_k$ after a time period $T$
  3. Using Arnoldi iteration calculate the new $\zeta_{k + 1}$, and complete the column k of the matrix $H$.
  
After $m$ iteration, it is possible calculate the $H$ eigenvalues. Finally the eigenvalues are scaled using the $\lambda_i = -log(\sigma_i)/T$. To reconstruct the eigenvector it is necessary to ise the columns of $H$ matrix.





