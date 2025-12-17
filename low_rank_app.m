
function runhist = low_rank_app (U,V,W)
    fprintf("low rank tensor approximation\n");
    timect=cputime;

    % the operator representaion of teosor A
    MA=khatrirao(U,V)*W';
    
    [m,r]=size(U);
    [n,r]=size(V);
    [p,r]=size(W);

    paras.m=m;
    paras.n=n;
    paras.p=p;
    paras.r=2;
    paras.itmax=150;
    rho=0.75;
    paras.tol=1e-4;

    gamma=1;

    runhist_vp = variable_projection_gamma(MA,rho,gamma,paras);
    %runhist_vp = variable_projection(MA, rho, paras,X);
    %runhist_vp = LM_method(MA,rho,gamma,paras);
    Uhat = runhist_vp.U;
    Vhat = runhist_vp.V;

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