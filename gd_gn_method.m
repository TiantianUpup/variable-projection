%%% this is a two-phase method
function runhist = gd_gn_method(U0,V0,MA,alpha,beta,epsilon)
    timect = cputime;
    [m,r]=size(U0);
    [n,r]=size(V0);
    U=U0;
    V=V0;
    runhist_gd=gradient_descent(U,V,MA,alpha,beta,epsilon);

    U=runhist_gd.U;
    V=runhist_gd.V;
    iter_gd=runhist_gd.iter;

    paras.m=m;
    paras.n=n;
    paras.r=r;
    paras.itmax=150;
    rho=0.75;
    paras.tol=1e-4;
    gamma=1;
    %%%%%%%%%%% variable projection method for gamma=1
    runhist_vp = variable_projection_gamma(U,V,MA,rho,gamma,paras);
    Uhat = runhist_vp.U;
    Vhat = runhist_vp.V;
    iter_gn=runhist_vp.iter;

    runhist.U=Uhat;
    runhist.V=Vhat;
    runhist.cput = cputime - timect;
    runhist.iter_gd=iter_gd;
    runhist.iter_gn=iter_gn;
end    