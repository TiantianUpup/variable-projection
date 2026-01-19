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
    cput_gd=runhist_gd.cput;

    paras.m=m;
    paras.n=n;
    paras.r=r;
    paras.itmax=150;
    paras.tol=1e-6;
    gamma=0;
    %%%%%%%%%%% variable projection method for gamma=1
    runhist_vp = variable_projection_gamma(U,V,MA,gamma,paras);
    Uhat = runhist_vp.U;
    Vhat = runhist_vp.V;
    iter_gn=runhist_vp.iter;
    cput_gn=runhist_vp.cput;


    runhist.U=Uhat;
    runhist.V=Vhat;
    runhist.cput = cputime - timect;
    runhist.iter_gd=iter_gd;
    runhist.iter_gn=iter_gn;
    runhist.cput_gn=cput_gn;
    runhist.cput_gd=cput_gd;
end    