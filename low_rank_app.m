
function runhist = low_rank_app (U,V,W)
    fprintf("low rank tensor approximation\n");
    timect=cputime;

    % the operator representaion of teosor A
    MA=khatrirao(U,V)*W';
    
    [m,r]=size(U);
    [n,r]=size(V);
    [p,r]=size(W);

    r=2;
    paras.m=m;
    paras.n=n;
    paras.r=r;
   

    vp_paras.gamma=0;
    vp_paras.itmax=1500;
    vp_paras.tol=1e-6; 
    vp_paras.mu=0;

    cg_paras.itmax=500;
    cg_paras.tol=1e-6;
    cg_paras.lambda=1;

    U = rand(m,r);
    V = rand(n,r);

    %runhist_vp = variable_projection_gamma(U,V,MA,vp_paras,cg_paras,paras);
    %runhist_vp = hybrid(U,V,MA,vp_paras,cg_paras,paras);
    runhist_vp = BFGS(U,V,MA,paras);
    %runhist_vp = LM_method(U,V,MA,vp_paras,cg_paras,paras);
    %runhist_vp = variable_projection(MA, rho, paras,X);
    %runhist_vp = LM_method(U,V,MA,gamma,paras);
    Uhat = runhist_vp.U;
    Vhat = runhist_vp.V;
    cput=runhist_vp.cput;
    fprintf("variable_projection method costs %3.4f\n",cput);

    kr=khatrirao(Uhat,Vhat);
    krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
    What=MA'*krmp';

    fprintf("====== Uhat ======\n");
    disp(Uhat);

    fprintf("====== Vhat ======\n");
    disp(Vhat);

    fprintf("====== What ======\n");
    disp(What);

    runhist.U=Uhat;
    runhist.V=Vhat;
    runhist.W=What;
    runhist.cput=cputime-timect;
end    