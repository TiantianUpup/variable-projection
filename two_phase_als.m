function runhist = two_phase_als(Uinit,X,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras)
    timect = cputime;
    % options_als = struct;
    % options_als.Compression = false;
    % options_als.Algorithm = @cpd_als;
    % options_als.AlgorithmOptions.MaxIter = 50;      % Default 500
    % options_als.AlgorithmOptions.CGMaxIter = 15;     % Default 15
    % options_als.Initialization = @cpd_rnd;
    % options_als.AlgorithmOptions.TolFun = 1e-3; 
    % options_als.AlgorithmOptions.TolX = 1e-3;
    
    r=paras.r;
    % [Xhat,out]=cpd(X,r,options_als);
    [Xhat,U0,out] = cp_als(tensor(X),r,...
        'printitn',0 , ...
        'init', Uinit, ...
        'maxiters', 150, ...
        'tol', 1e-3);

    runhist.cput_als = cputime-timect;
    runhist.iter_als = out.iters;

    fprintf("---------- two_phase_als method, als iter=%d ---------\n", out.iters);

    U0=Xhat{1};
    V0=Xhat{2};

    % variable projection method
    runhist_vp = variable_projection_simple(U0,V0,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);

    runhist.U = runhist_vp.U;
    runhist.V = runhist_vp.V;
    runhist.W = runhist_vp.W;
    runhist.iter_vp = runhist_vp.iter;
    runhist.cput_vp = runhist_vp.cput;
    runhist.cput = cputime - timect;
    runhist.iter=runhist_vp.iter;
    runhist.fval=runhist_vp.fval;
end     