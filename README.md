# A Projected Gauss-Newton Variable Projection Method for Low Rank Approximations of Tensors



### Copyright
### Quickstart
`vp_pGN` is the core solver in the folder `algo`.

#### usage
```matlab
runhist = vp_pGN(T, Uinit, vp_paras, cg_paras);
```
- `T`: Factor matrices of the ground-truth tensor.
- `Uinit`: Initial factor matrices.
- `vp_paras`: Parameters of the variable projection method.
- `cg_paras`: Parameters of the conjugate gradient method.

The solver returns `runhist` containing:
- `runhist.U`, `runhist.V`, `runhist.W`: Computed factor matrices.
- `runhist.fval`: Objective function values during iterations.
- `runhist.iter`: Number of iterations.
- `runhist.cput`: Computational time.

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

All the numerical results are saved in the folder `results`

If you have any questions about implementation, please contact Tiantian He (hetiantian@nudt.edu.cn) and  Shenglong Hu(hushenglong@nudt.edu.cn)

