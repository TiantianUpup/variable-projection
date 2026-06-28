%%% this is a two-phase method
function runhist = gd_gn_method(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras)
    fprintf("==================== This is gd_gn_method ======================\n");
    timect = cputime;
    [m,r]=size(U);
    [n,r]=size(V);
   

    alpha=1e-4;
    beta=0.8;
    epsilon=1e-2;
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
    runhist_vp = variable_projection_simple(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    Uhat = runhist_vp.U;
    Vhat = runhist_vp.V;
    What = runhist_vp.W;
    iter_gn=runhist_vp.iter;
    cput_gn=runhist_vp.cput;

    runhist.U=Uhat;
    runhist.V=Vhat;
    runhist.W=What;
    runhist.cput = cputime - timect;
    runhist.iter_gd=iter_gd;
    runhist.iter=iter_gn;
    runhist.iter_gn=iter_gn;
    runhist.cput_gn=cput_gn;
    runhist.cput_gd=cput_gd;
    runhist.fval=runhist_vp.fval;
end    