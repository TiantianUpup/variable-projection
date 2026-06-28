clc; clear all; close all;


addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
addpath("E:\\matlab-code\\variable-projection\\tensorlab")
% m=30;
% n=40;
% p=1000;

% m=15;
% n=15;
% p=15;

% m=30;
% n=40;
% p=50;

% m=30;
% n=30;
% p=40;

m=100;
n=100;
p=100;

% m=30;
% n=40;
% p=1000;

% m=100;
% n=100;
% p=10000;

% m=30;
% n=30;
% p=40;

% m=2; 
% n=3;
% p=3;
trial=10;

%%%%%%%% variable projection method
iter_vp_result=zeros(trial,1);
cput_vp_result=zeros(trial,1);
res_vp_result=zeros(trial,1);
rel_error_vp_result=zeros(trial,1);

%%%%%%%% variable projection method (random initialization)
iter_vp_rand_result=zeros(trial,1);
cput_vp_rand_result=zeros(trial,1);
res_vp_rand_result=zeros(trial,1);
rel_error_rand_vp_result=zeros(trial,1);

%%%%%%%% v4
iter_vp_rand_result_4=zeros(trial,1);
cput_vp_rand_result_4=zeros(trial,1);
res_vp_rand_result_4=zeros(trial,1);
rel_error_rand_vp_result_4=zeros(trial,1);

%%%%%%%% variable projection gamma method
iter_vp_gamma_result=zeros(trial,1);
cput_vp_gamma_result=zeros(trial,1);
res_vp_gamma_result=zeros(trial,1);
rel_error_vp_gamma_result=zeros(trial,1);

%%%%%%%% cp_als method
iter_cp_als_result=zeros(trial,1);
cput_cp_als_result=zeros(trial,1);
res_cp_als_result=zeros(trial,1);
rel_error_cp_als_result=zeros(trial,1);

%%%%%%%% cpd_als method
iter_cpd_als_result=zeros(trial,1);
cput_cpd_als_result=zeros(trial,1);
res_cpd_als_result=zeros(trial,1);
rel_error_cpd_als_result=zeros(trial,1);

%%%%%%%% cpd_nls method
iter_cpd_nls_result=zeros(trial,1);
cput_cpd_nls_result=zeros(trial,1);
res_cpd_nls_result=zeros(trial,1);
rel_error_cpd_nls_result=zeros(trial,1);

c=0.9;

% U_true=rand(m,r);
% V_true=rand(n,r);
% W_true=rand(p,r);
Rel_err=zeros(10,4);
T=zeros(10,4);
count=1;
for i=1:trial
    % truth factor matrices
    r=15;
    U_true=rand(m,r);
    V_true=rand(n,r);
    W_true=rand(p,r);

    MA=khatrirao(U_true,V_true)*W_true';
    X=generate_cp_tensor(U_true,V_true,W_true);

    % N=randn(m,n,p);
    % rho_t=0.1*(norm(X(:))/norm(N(:)));
    % X_noise=X+rho_t*N;/norm(N(:))));
    % MA_noise=MA+rho*N;

%     kappa=1e6;
%     Z = orth(randn(p, r));
%     R = orth(randn(r, r));
%     s = logspace(0, -log10(kappa), r);
%     W_true = Z*diag(s)*R';
    % X=generate_cp_tensor(U_true,V_true,W_true);

    % % colinear factor matrices
    % T=col_mat_gen(m,r,c);
    % U_true=T{1};
    % V_true=T{2};
    % W_true=T{3};
    % X=generate_cp_tensor(U_true,V_true,W_true);
    % MA=khatrirao(U_true,V_true)*W_true';

    % %% complete orthogonality
    % U_true = orth(rand(m,r));
    % V_true = orth(rand(n,r));
    % W_true = orth(rand(p,r));


    %%%%%%%%%%% test data
    r=16;
    paras.m=m; 
    paras.n=n;   
    paras.r=r;  

    % U = orth(rand(m,r));
    % V = orth(rand(n,r));

    % initial point
    U = rand(m,r);
    U = proj_oblique(U);
    V = rand(n,r);
    V = proj_oblique(V);
    W = rand(p,r);
    W = proj_oblique(W);
    Uinit={U,V,W};

    % gevd initialization
    % cput_gevd=cputime;
    % [Uinit,output] = cpd_gevd(X,r);
    % fprintf("--------------------------- cpd_gevd costs %3.10f ----------------------------\n",cputime-cput_gevd);

    % U=Uinit{1};
    % V=Uinit{2};
    
    % fprintf("the size of U is\n");
    % disp(size(U));

    vp_paras.gamma=0;
    vp_paras.itmax=1500;
    vp_paras.tol=1e-6; 
    vp_paras.mu=0;
    
    cg_paras.itmax=500;
    cg_paras.tol=1e-6;
    cg_paras.lambda=0;

    %profile on;
    runhist_vp_0 = variable_projection_simple_opt_v5(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    %runhist_vp_0 = variable_projection_simple(U,V,U_true,V_true,W_true,MA_noise,vp_paras,cg_paras,paras);
    %profile off;
    %profile viewer;
    Uhat_1 = runhist_vp_0.U;
    Vhat_1 = runhist_vp_0.V;
    What_1 = runhist_vp_0.W;
    fval_vp = runhist_vp_0.fval;
    % cput_vp_0 = runhist_vp_0.cput;
    cput_vp_0 = runhist_vp_0.cput;
    iter_vp_0 = runhist_vp_0.iter;

    Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
    Xres0=Xhat_1-X;

    iter_vp_result(i)=iter_vp_0;
    cput_vp_result(i)=cput_vp_0;
    res_vp_result(i)=norm(Xres0(:));
    %rel_error_vp_result(i)=norm(Xres0(:))/norm(X(:));
    rel_error_vp_result(i)=norm(Xres0(:))/norm(X(:));

    % % % gevd initialization
    % % cput_gevd=cputime;
    % % [Uinit_gevd,output] = cpd_gevd(X,r);
    % % fprintf("--------------------------- cpd_gevd costs %3.10f ----------------------------\n",cputime-cput_gevd);

    % % U_rand=Uinit_gevd{1};
    % % V_rand=Uinit_gevd{2};
    
    % runhist_vp_rand = variable_projection_simple_opt_v4(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % %profile off;
    % %profile viewer;
    % Uhat_vp_rand = runhist_vp_rand.U;
    % Vhat_vp_rand = runhist_vp_rand.V;
    % What_vp_rand = runhist_vp_rand.W;
    % fval_vp_rand = runhist_vp_rand.fval;
    % cput_vp_rand = runhist_vp_rand.cput;
    % iter_vp_rand = runhist_vp_rand.iter;

    % Xhat_rand=generate_cp_tensor(Uhat_vp_rand,Vhat_vp_rand,What_vp_rand);
    % Xres_vp_rand=Xhat_rand-X;

    % iter_vp_rand_result(i)=iter_vp_rand;
    % cput_vp_rand_result(i)=cput_vp_rand;
    % res_vp_rand_result(i)=norm(Xres_vp_rand(:));
    % rel_error_rand_vp_result(i)=norm(Xres_vp_rand(:))/norm(X(:));

    % % runhist_vp_rand = two_phase_als(Uinit,X,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);

    % %runhist_vp_rand = variable_projection_simple_opt(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % runhist_vp_rand = variable_projection_simple_opt_v3(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % %runhist_vp_rand = variable_projection_simple(U,V,U_true,V_true,W_true,MA_noise,vp_paras,cg_paras,paras);
    % %profile off;
    % %profile viewer;
    % Uhat_vp_rand = runhist_vp_rand.U;
    % Vhat_vp_rand = runhist_vp_rand.V;
    % What_vp_rand = runhist_vp_rand.W;
    % fval_vp_rand = runhist_vp_rand.fval;
    % cput_vp_rand = runhist_vp_rand.cput;
    % iter_vp_rand = runhist_vp_rand.iter;

    % Xhat_rand=generate_cp_tensor(Uhat_vp_rand,Vhat_vp_rand,What_vp_rand);
    % Xres_vp_rand=Xhat_rand-X;

    % iter_vp_rand_result(i)=iter_vp_rand;
    % cput_vp_rand_result(i)=cput_vp_rand;
    % res_vp_rand_result(i)=norm(Xres_vp_rand(:));
    % rel_error_rand_vp_result(i)=norm(Xres_vp_rand(:))/norm(X(:));

    % %%%%%%%%%%%%%%%%%%%%% version 4 code
    % runhist_vp_rand_4 = variable_projection_simple_opt_v4(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % %runhist_vp_rand = variable_projection_simple(U,V,U_true,V_true,W_true,MA_noise,vp_paras,cg_paras,paras);
    % %profile off;
    % %profile viewer;
    % Uhat_vp_rand_4 = runhist_vp_rand_4.U;
    % Vhat_vp_rand_4 = runhist_vp_rand_4.V;
    % What_vp_rand_4 = runhist_vp_rand_4.W;
    % fval_vp_rand_4 = runhist_vp_rand_4.fval;
    % cput_vp_rand_4 = runhist_vp_rand_4.cput;
    % iter_vp_rand_4 = runhist_vp_rand_4.iter;

    % Xhat_rand_4=generate_cp_tensor(Uhat_vp_rand_4,Vhat_vp_rand_4,What_vp_rand_4);
    % Xres_vp_rand_4=Xhat_rand_4-X;

    % iter_vp_rand_result_4(i)=iter_vp_rand_4;
    % cput_vp_rand_result_4(i)=cput_vp_rand_4;
    % res_vp_rand_result_4(i)=norm(Xres_vp_rand_4(:));
    % rel_error_rand_vp_result_4(i)=norm(Xres_vp_rand_4(:))/norm(X(:));



    
    % vp_paras.gamma=1;
    
    % runhist_vp_0 = variable_projection_simple(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % %profile off;
    % %profile viewer;
    % Uhat_1_gamma = runhist_vp_0.U;
    % Vhat_1_gamma = runhist_vp_0.V;
    % What_1_gamma = runhist_vp_0.W;
    % fval_vp_gamma =runhist_vp_0.fval;
    % cput_vp_gamma = runhist_vp_0.cput;
    % iter_vp_gamma = runhist_vp_0.iter;

    % Xhat_1=generate_cp_tensor(Uhat_1_gamma,Vhat_1_gamma,What_1_gamma);
    % Xres_gamma=Xhat_1-X;

    % iter_vp_gamma_result(i)=iter_vp_gamma;
    % cput_vp_gamma_result(i)=cput_vp_gamma;
    % res_vp_gamma_result(i)=norm(Xres_gamma(:));
    % rel_error_vp_gamma_result(i)=norm(Xres_gamma(:))/norm(X(:));


    % % runhist_mcg = variable_projection_simple_acg(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % runhist_mcg = variable_projection_simple_acg(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
    % Uhat_1 = runhist_mcg.U;
    % Vhat_1 = runhist_mcg.V;
    % What_1 = runhist_mcg.W;
    % cput_mcg = runhist_mcg.cput;
    % iter_mcg = runhist_mcg.iter;

    % Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
    % Xres_mcg=Xhat_1-X;

    % tensor_toolbox method
    timect = cputime;
    [P,U0,out] = cp_als(tensor(X),r,...
        'printitn',0 , ...
        'init', Uinit, ...
        'maxiters', 1500, ...
        'tol', 1e-6);

    % cput_als = cputime - timect;   
    cput_als = cputime - timect;   
    als_res=full(P)-X;
    fval_als=out.fval;
    iter_als=out.iters;

    iter_cp_als_result(i)=iter_als;
    cput_cp_als_result(i)=cput_als;
    res_cp_als_result(i)=norm(als_res(:));
    %rel_error_cp_als_result(i)=norm(als_res(:))/norm(X(:));
    rel_error_cp_als_result(i)=norm(als_res(:))/norm(X(:));

    % tensorlab method
    options_als = struct;
    %options_als.Compression = false;
    options_als.Algorithm = @cpd_als;
    options_als.AlgorithmOptions.MaxIter = 1500;      % Default 500
    timect = cputime;
    [Xhat,out_cpd_als]=cpd(X,Uinit,options_als);
    %[Xhat,out_cpd_als]=cpd(X_noise,Uinit,options_als);
    %[Xhat,out_cpd_als]=cpd(X,r,options_als);
    %cput_cpd_als = cputime - timect;
    cput_cpd_als = cputime - timect;
    lab_res=X-cpdgen(Xhat);
    fval_cpd_als=out_cpd_als.Algorithm.fval;
    iter_cpd_als=out_cpd_als.Algorithm.iterations;
    
    disp(out_cpd_als.Algorithm);
    disp(out_cpd_als.Compression);
    disp(out_cpd_als.Preprocessing);

    iter_cpd_als_result(i)=iter_cpd_als;
    cput_cpd_als_result(i)=cput_cpd_als;
    res_cpd_als_result(i)=norm(lab_res(:));
    % rel_error_cpd_als_result(i)=norm(lab_res(:))/norm(X(:));
    rel_error_cpd_als_result(i)=norm(lab_res(:))/norm(X(:));

    % nls method
    options_nls = struct;
    options_nls.Compression = false;
    options_nls.Algorithm = @cpd_nls;
    options_nls.AlgorithmOptions.MaxIter = 1500;      % Default 500
    options_nls.AlgorithmOptions.CGMaxIter = 500;      % Default 15
    timect_rnd = cputime;
    
    [Uhat,out_cpd_nls]=cpd(X,Uinit,options_nls);
    %[Uhat,out_cpd_nls]=cpd(X_noise,Uinit,options_nls);
    %[Uhat,out_cpd_nls]=cpd(X,r,options_nls);
    %cput_cpd_nls = cputime - timect_rnd;
    cput_cpd_nls = cputime - timect_rnd;
    lab_res_rand=X-cpdgen(Uhat);
    fval_cpd_nls=out_cpd_nls.Algorithm.fval;
    iter_cpd_nls=out_cpd_nls.Algorithm.iterations;

    disp(out_cpd_nls.Algorithm);
    disp(out_cpd_nls.Compression);
    disp(out_cpd_nls.Preprocessing);

    iter_cpd_nls_result(i)=iter_cpd_nls;
    cput_cpd_nls_result(i)=cput_cpd_nls;
    res_cpd_nls_result(i)=norm(lab_res_rand(:));
    % rel_error_cpd_nls_result(i)=norm(lab_res_rand(:))/norm(X(:));
    rel_error_cpd_nls_result(i)=norm(lab_res_rand(:))/norm(X(:));

    Rel_err(count,:)=[rel_error_cp_als_result(i),rel_error_cpd_als_result(i),rel_error_cpd_nls_result(i),rel_error_vp_result(i)];
    T(count,:)=[cput_als,cput_cpd_als,cput_cpd_nls,cput_vp_0];
    count=count+1;

    fprintf("============================= iter=%d ===================\n",i);
    fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f, res_err is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)),norm(Xres0(:))/norm(X(:)));
    %fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (random) is %3.10f, res_err is %3.10f\n",iter_vp_rand,cput_vp_rand,norm(Xres_vp_rand(:)),norm(Xres_vp_rand(:))/norm(X(:)));
    %fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (random) is %3.10f, res_err is %3.10f\n",iter_vp_rand_4,cput_vp_rand_4,norm(Xres_vp_rand_4(:)),norm(Xres_vp_rand_4(:))/norm(X(:)));
    %fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f, res_err is %3.10f\n",iter_vp_gamma,cput_vp_gamma,norm(Xres_gamma(:)),norm(Xres_gamma(:))/norm(X(:)));
    %fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_mcg (gamma=0) is %3.10f, res_err is %3.10f\n",iter_mcg,cput_mcg,norm(Xres_mcg(:)),norm(Xres_mcg(:))/norm(X(:)));
    fprintf("iter=%d, cput=%3.4f, the residual of the als method is %3.10f, res_err is %3.10f\n",iter_als,cput_als,norm(als_res(:)),norm(als_res(:))/norm(X(:)));
    fprintf("iter=%d, cput=%3.4f, the residual of the cpd_als method is %3.10f, res_err is %3.10f\n",iter_cpd_als,cput_cpd_als,norm(lab_res(:)), norm(lab_res(:))/norm(X(:)));
    fprintf("iter=%d, cput=%3.4f, the residual of the cpd_nls method is %3.10f, res_err is %3.10f\n",iter_cpd_nls,cput_cpd_nls,norm(lab_res_rand(:)),norm(lab_res_rand(:))/norm(X(:)));
end

fprintf("========================================== rank %d approximation ==========================================\n",r);

%fprintf("%d &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f\n",r, mean(iter_vp_result),mod(mean(res_vp_result),1),mean(cput_vp_result),mean(iter_cp_als_result),mod(mean(res_cp_als_result),1),mean(cput_cp_als_result),mean(iter_cpd_als_result),mod(mean(res_cpd_als_result),1),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mod(mean(res_cpd_nls_result),1),mean(cput_cpd_nls_result));
fprintf("%d &%3.1f &%.4f &%3.1f &%3.1f &%.4f &%3.1f &%3.1f &%.4f &%3.1f &%3.1f &%.4f &%3.1f\n",r,mean(iter_cp_als_result),mean(res_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(res_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(res_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(res_vp_result),mean(cput_vp_result));
%fprintf("%d &%3.1f &%.4f &%3.1f &%3.1f &%.4f &%3.1f &%3.1f &%.4f &%3.1f &%3.1f &%.4f &%3.1f\n",r,mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result));
fprintf("%d &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f\n",r,mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result));

% fprintf("relative error of vp=%.4f, vp-gamma=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(rel_error_vp_result),mean(rel_error_vp_gamma_result),mean(rel_error_cp_als_result),mean(rel_error_cpd_als_result),mean(rel_error_cpd_nls_result));
% fprintf("mean value of residul value of vp=%.4f,vp-gamma=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(res_vp_result),mean(res_vp_gamma_result),mean(res_cp_als_result),mean(res_cpd_als_result),mean(res_cpd_nls_result));

% fprintf("relative error of vp=%.4f, vp_rand=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(rel_error_vp_result),mean(rel_error_rand_vp_result),mean(rel_error_cp_als_result),mean(rel_error_cpd_als_result),mean(rel_error_cpd_nls_result));
% fprintf("mean value of residul value of vp=%.4f, vp_rand=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(res_vp_result),mean(res_vp_rand_result),mean(res_cp_als_result),mean(res_cpd_als_result),mean(res_cpd_nls_result));

fprintf("relative error of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(rel_error_cp_als_result),mean(rel_error_cpd_als_result),mean(rel_error_cpd_nls_result),mean(rel_error_vp_result));
fprintf("mean value of residul value of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(res_cp_als_result),mean(res_cpd_als_result),mean(res_cpd_nls_result),mean(res_vp_result));

% % FIVE test case
% fprintf("mean value of iteration of vp=%.4f, vp_rand=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(iter_vp_result),mean(iter_vp_rand_result),mean(iter_cp_als_result),mean(iter_cpd_als_result),mean(iter_cpd_nls_result));
% fprintf("mean value of cpu time of vp=%.4f, vp_rand=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(cput_vp_result),mean(cput_vp_rand_result),mean(cput_cp_als_result),mean(cput_cpd_als_result),mean(cput_cpd_nls_result));

% fprintf("mean value of iteration of vp=%.4f, vp_rand=%.4f, vp_rand_4=%.4f,cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(iter_vp_result),mean(iter_vp_rand_result),mean(iter_vp_rand_result_4),mean(iter_cp_als_result),mean(iter_cpd_als_result),mean(iter_cpd_nls_result));
% fprintf("mean value of cpu time of vp=%.4f, vp_rand=%.4f, vp_rand_4=%.4f, cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f\n",mean(cput_vp_result),mean(cput_vp_rand_result),mean(cput_vp_rand_result_4),mean(cput_cp_als_result),mean(cput_cpd_als_result),mean(cput_cpd_nls_result));

fprintf("mean value of iteration of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(iter_cp_als_result),mean(iter_cpd_als_result),mean(iter_cpd_nls_result),mean(iter_vp_result));
fprintf("mean value of cpu time of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(cput_cp_als_result),mean(cput_cpd_als_result),mean(cput_cpd_nls_result),mean(cput_vp_result));

% iter_vp_rand_result(i)=iter_vp_rand_4;
% cput_vp_rand_result(i)=cput_vp_rand_4;


fprintf("vp, min=%3.4f, max=%3.4f\n",min(rel_error_vp_result),max(rel_error_vp_result));
fprintf("cp_als, min=%3.4f, max=%3.4f\n",min(rel_error_cp_als_result),max(rel_error_cp_als_result));
fprintf("cpd_als, min=%3.4f, max=%3.4f\n",min(rel_error_cpd_als_result),max(rel_error_cpd_als_result));
fprintf("cpd_nls, min=%3.4f, max=%3.4f\n",min(rel_error_cpd_nls_result),max(rel_error_cpd_nls_result));

Result = [res_cp_als_result-res_vp_result, res_cpd_als_result-res_vp_result, res_cpd_nls_result-res_vp_result];
prob=sum(all(Result >= 0, 2));
% prob=
max_gap=max(max(Result(Result >= 0)));
fprintf("the probability is %3.1f, the max residual is %3.2f\n",prob/trial,max_gap);
disp("Result matrix is:");
disp(Result);

disp("Rel_err matrix is:");
disp(Rel_err);

formatSpec = '%.4f %.4f %.4f %.4f\n';
r_id = fopen('rel_err.txt', 'w');
t_id = fopen('time.txt', 'w');
fprintf(r_id, formatSpec, Rel_err');
fprintf(t_id, formatSpec, T');

fclose(r_id); 
fclose(t_id); 
disp("Time matrix is:");
disp(T);

% fprintf(r, mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result),mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result));

% Result_res = [res_cp_als_result-res_vp_result, res_cpd_als_result-res_vp_result, res_cpd_nls_result-res_vp_result];
% prob_res=sum(all(Result_res >= 0, 2));
% % prob=
% max_gap=max(max(Result_res(Result_res >= 0)));
% fprintf("the probability is %3.1f, the max residual is %3.2f\n",prob_res/trial,max_gap);
% disp("Result_res matrix is:");
% disp(Result_res);

fprintf("rho is %d\n",rho); 

 