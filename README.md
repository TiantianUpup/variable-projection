# A Projected Gauss-Newton Variable Projection Method for Low Rank Approximations of Tensors
This package includes a MATLAB implementation of the algorithm vp_pGN presented in the paper "A projected Gauss-Newton variable projection method for low rank approximations of tensors" by Tiantian He, Shenglong Hu and Zheng-Hai Huang.

### Introduction
The `algo` folder contains our core solver `vp_pGN`. The solver requires the following parameters:
```matlab
runhist = vp_pGN(T, Uinit, vp_paras, cg_paras);
```
- `T`: factor matrices of the ground-truth tensor.
- `Uinit`: initial factor matrices.
- `vp_paras`: parameters of the variable projection method.
- `cg_paras`: parameters of the conjugate gradient method.

The solver returns `runhist` containing:
- `runhist.U`, `runhist.V`, `runhist.W`: computed factor matrices.
- `runhist.fval`: the history of the relative changes of the objective function.
- `runhist.iter`: number of iterations.
- `runhist.cput`: CPU time.

For example:
```matlab
m=50; n=50; p=50;
r_true=15;

U_true=rand(m,r_true); V_true=rand(n,r_true);
W_true=rand(p,r_true); 

T={U_true,V_true,W_true};
X=generate_cp_tensor(U_true,V_true,W_true);

r=10; % approximation rank

U = rand(m,r); U = proj_oblique(U);
V = rand(n,r); V = proj_oblique(V);
W = rand(p,r); W = proj_oblique(W);
Uinit={U,V,W};

vp_paras.itmax=1500; vp_paras.xtol=1e-6; vp_paras.ftol=1e-12; 
cg_paras.itmax=500; cg_paras.tol=1e-6;

runhist = vp_pGN(T, Uinit, vp_paras, cg_paras);
```

### Quickstart
The folder `examples` contains five presented in the paper.
- `test_exm1.m`: reproduces the numerical experiments in Subsection 6.1 and tests three problem sizes: $m=n=p=100$, $m=30,n=40, p=1000$, and $m=n=100,p=10000$.
- `test_exm2.m`: 
- `test_exm3.m`:
- `test_exm4.m`:
- `test_exm5.m`:

All the numerical results are saved in the folder `results`.

If you have any questions about implementation, please contact Tiantian He (hetiantian@nudt.edu.cn) and  Shenglong Hu(hushenglong@nudt.edu.cn).

