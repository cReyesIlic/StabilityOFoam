# StabilityOFoam 1.0

## Solvers 
Solvers and cases to run a Jacobian-free matrix linear stability analysis.
Tested in OpenFOAM 7.x with cylinder wake (stable and unstable cases)

Steps for compiling
1. Install OpenFOAM (visit https://openfoam.org/)
2. Compiling icoFOAMSFD solver
3. Compiling icoFOAMLNS solver

The OF implementation is based on paper found on https://www.sciencedirect.com/science/article/pii/S1000936116300024.

### icoFOAMSFD solver

Solver icoFOAM based to calculate the base flow for stability analysis using selective frequency damping (SFD). It is necessary to define in transportProperties dictionary the valu of $\chi$ and rDelta (parameter of SFD).

#### Case Re = 100

Complete damping: $\chi = 0.5$, $1/\Delta = 0.5$

![alt text](https://github.com/cReyesIlic/StabilityOFoam/blob/main/img/Xi_rDelta05.png?raw=true)

Partial damping: $\chi = 2$, $1/\Delta = 0.5$

![alt text](https://github.com/cReyesIlic/StabilityOFoam/blob/main/img/rDelta2Xi05.png?raw=true)

### icoFOAMLNS

solver icoFOAM based which calculate the perturbation velocity by Linearized Navier Stokes,

## Time-stepper approach

The system linearizated can be expressed as

$\frac{\partial u'}{t}  =Au'$

where $A$ is the projected Jacobian operator. A solution for this problem corresponds to $u'(\Delta t) = u'_0 e^{A\Delta t} = u'(\Delta t)B$, where $e^{A\delta t}$ is called the exponential propagator where matrix A eigenpairs is related to B matrix eigenpairs by $\Lambda = \log(\Sigma)/\Delta t$ and $V = V_e$, where $V_e$ is obtained from the Hessenberg Matrix obtained by a Arnoldi iteration.

## Linear Stability Analysis - Jacobian Free Method- Hessenberg Matrix

To calculate the H matrix (size m x m, with m and user-specified parameter) it is necessary to follow the next steps:

1. Calculate the initial random $\zeta$ value, where $\zeta = [U_1\quad V_1 \quad W_1 \quad U_2 \quad V_2 \quad W_2 \quad ...]^T $
Begin Iteration k = 0:
  2. Calculate the new perturbation $U_{k + 1} = \beta\zeta_k$ after a time period $T$
  3. Using Arnoldi iteration calculate the new $\zeta_{k + 1}$, and complete the column k of the matrix $H$.
  
After $m$ iteration, it is possible calculate the $H$ eigenvalues. Finally the eigenvalues are scaled using the $\lambda_i = -log(\sigma_i)/T$. To reconstruct the eigenvector it is necessary to ise the columns of $H$ matrix.

## Example 1 : Cylinder wake Re = 40 (Stable case)

https://github.com/cReyesIlic/StabilityOFoam/blob/main/img/Re40.png
![alt text](https://github.com/cReyesIlic/StabilityOFoam/blob/main/img/Re40.png?raw=true)

## Example 2: Cylinder wake Re = 100 (Unstable case)
![alt text](https://github.com/cReyesIlic/StabilityOFoam/blob/main/img/Re100.png?raw=true)




