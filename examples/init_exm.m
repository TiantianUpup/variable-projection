clc; clear all; close all;

addpath('../compared_method/tensor_toolbox');
addpath('../compared_method/tensorlab');
addpath('../utils');
addpath('../algo');


trial=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Running history
%%%%%%%% vp_pGN method
iter_vp_result=zeros(trial,1);
cput_vp_result=zeros(trial,1);
res_vp_result=zeros(trial,1);
rel_error_vp_result=zeros(trial,1);

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

Rel_err=zeros(10,4);
Time=zeros(10,4);
count=1;

m=50;
n=50;
p=50;

kappa=1e4;
%for noise=[5,10,15,20,25,30,35,40,45,55,60,65,70,75,80]
%for epsion=[0.1,0.5,0.01,0.05,0.001,0.005]
epsilon=1e-2;
for r=1:15
for i=1:trial
    r_true=15;
    % [R,~] = qr(rand(r_true),0);  
    % [Q1,~] = qr(rand(m,r_true),0);
    % [Q2,~] = qr(rand(n,r_true),0);
    % [Q3,~] = qr(rand(p,r_true),0);
    % sigma = logspace(0,-log10(kappa),r_true);
    % Sigma = diag(sigma);

    % sigma_1 = logspace(0,-log10(kappa),r_true);
    % Sigma_1 = diag(sigma_1);

    % U_true = Q1*Sigma*R';
    % %V_true = Q2*Sigma*R';
    % V_true=rand(n,r_true);
    % W_true = Q3*Sigma_1*R';


    % N = randn(size(V_true));
    % DeltaA = epsilon * norm(V_true(:)) ...
    %    / norm(N(:)) * N;
    % V_noise=V_true+DeltaA;
    % C = 3;
    % eta = 0.5;
    % % lambda = zeros(n,1);
    % X=generate_cp_tensor(U_true,V_true,W_true);
    % for kk = 1:r_true
    %     lam=C*norm(X(:),'fro')/(kk^(1+eta));
    %     W_temp(:,kk)=lam*W_temp(:,kk);
    % end

    
    %Lambda = diag(lambda);

    % W_noise=W_true+0.01*W_temp;
  
    % [Q,~] = qr(randn(p,r_true),0);
    % [R,~] = qr(randn(r_true,r_true));
    % eps = 1e-4;
    % s = [ones(r_true-1,1); eps];
    % W_true = Q*diag(s)*R';


    % [U_true,V_true,W_true]=generate_factors(m, n, p, r_true);

    % % AA = rand(p,14);
    % % BB = rand(r_true,14);
    % % W_true = AA*BB';
    % W_true=rand(p,r_true);
    % %fprintf("==================== the rank of the W_true is %d ===============\n",rank(W_true));
    % noise_level = 1e-3;
    % W_noise = W_true + noise_level*rand(p,r_true);
    % %T={U_true,V_true,W_true};
    
    % density = 0.1;
    % W_true=rand(p,r_true);
    % S = zeros(size(W_true));
    % idx = randperm(numel(W_true), ...
    %     round(density*numel(W_true)));
    % S(idx)=randn(length(idx),1);
    % rho = 0.1; % relative perturbation
    % DeltaW = rho * norm(W_true,'fro') ...
    %         /norm(S,'fro') * S;
    % W_noise = W_true + DeltaW;

    % delta=0.001;
    % W_true=rand(p,r_true);
    % % N = randn(size(W_true));
    % % E = delta * N / norm(N,2);
    % % W_noise = W_true + E;
    
    % [W_noise,N]=noisy(W_true,noise);

    
    % W_n=randn(p,r_true);
    % rho_w=epsion*(norm(W_true(:))/norm(W_n(:)));
    % % % disp(rho_w);
    % W_noise=W_true+rho_w*W_n;
    

    %[W_noise,N]=noisy(W_true,noise);

    %fprintf("==================== the rank of the W_noise is %d ===============\n",rank(W_noise));


    U_true=rand(m,r_true);
    V_true=rand(n,r_true);
    W_true=rand(p,r_true);
    % lam=[15 10 7 5 3 ...
    %       2 1.3 0.8 0.5 0.3 ...
    %       0.2 0.12 0.08 0.05 0.03];
    % lam=[80 50 30 18 10 ...
    %       5 2.5 1.2 0.6 0.3 ...
    %       0.15 0.08 0.04 0.002 0.001];      
    
    lam=logspace(1,-1,15);

    %lam=[100,50,20,10,5,2,1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.001];

    W_true=W_true*diag(lam);
    % kappa=1e6;
    % Z = orth(randn(p, r_true));
    % R = orth(randn(r_true, r_true));
    % s = logspace(0, -log10(kappa), r_true);
    % W_true = Z*diag(s)*R';

    % W_n=randn(p,r_true);
    % % 0.001
    % rho_w=0.001*(norm(W_true(:))/norm(W_n(:)));
    % % % disp(rho_w);
    % W_noise=W_true+rho_w*W_n;

    % [Q,R] = qr(W_true,0);
    % epsilon = 1e-3;
    % D = epsilon*diag(randn(r_true,1));
    % R_new = R + D;
    % W_noise = Q*R_new;

    %T={U_true,V_true,W_noise};
    T={U_true,V_true,W_true};
    X=generate_cp_tensor(U_true,V_true,W_true);
    %X_noise=generate_cp_tensor(U_true,V_true,W_noise);
    X_noise=generate_cp_tensor(U_true,V_true,W_true);

    % %%%%%%%%%%% test data
    % r=12;
    % paras.m=m; 
    % paras.n=n;   
    % paras.r=r; 

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
    Uinit_vp={U,V};

    % gevd initialization
    % cput_gevd=cputime;
    % [Uinit,output] = cpd_gevd(X,r);
    % fprintf("--------------------------- cpd_gevd costs %3.10f ----------------------------\n",cputime-cput_gevd);

    % U=Uinit{1};
    % V=Uinit{2};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% vp_pGN method (ours)
    vp_paras.itmax=1500;
    vp_paras.xtol=1e-6; 
    vp_paras.ftol=1e-12; 
    
    cg_paras.itmax=500;
    cg_paras.tol=1e-6;

    % options_vp = struct;
    % options_vp.NLS.xtol = 1e-6; 
    % options_vp.NLS.ftol = 1e-12;
    % options_vp.NLS.itmax = 1500;

    % options_vp.CG.itmax = 500;
    % options_vp.CG.tol = 1e-6;

    runhist_vp_0 = vp_pGN(T, Uinit, vp_paras, cg_paras);
    Uhat_1 = runhist_vp_0.U;
    Vhat_1 = runhist_vp_0.V;
    What_1 = runhist_vp_0.W;
    fval_vp = runhist_vp_0.fval;
    cput_vp_0 = runhist_vp_0.cput;
    iter_vp_0 = runhist_vp_0.iter;

    Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
    Xres0=Xhat_1-X;

    iter_vp_result(i)=iter_vp_0;
    cput_vp_result(i)=cput_vp_0;
    res_vp_result(i)=norm(Xres0(:));
    rel_error_vp_result(i)=norm(Xres0(:))/norm(X(:));


    %%%%%%%%%%%%%%%%%%%%%%%%% cp_als method
    timect = cputime;
    [P,U0,out] = cp_als(tensor(X_noise),r,...
        'printitn',0 , ...
        'init', Uinit, ...
        'maxiters', 1500, ...
        'tol', 1e-8);

    cput_als = cputime - timect;   
    als_res=full(P)-X;
    fval_als=out.fval;
    iter_als=out.iters;

    iter_cp_als_result(i)=iter_als;
    cput_cp_als_result(i)=cput_als;
    res_cp_als_result(i)=norm(als_res(:));
    rel_error_cp_als_result(i)=norm(als_res(:))/norm(X(:));

    
    %%%%%%%%%%%%%%%%%%%%%%%%% cpd_als method
    options_als = struct;
    options_als.Algorithm = @cpd_als;
    options_als.AlgorithmOptions.MaxIter = 1500;      
    options_als.AlgorithmOptions.TolFun = 1e-12; 
    timect = cputime;
    [Xhat,out_cpd_als]=cpd(X_noise,Uinit,options_als);
    cput_cpd_als = cputime - timect;
    lab_res=X-cpdgen(Xhat);
    fval_cpd_als=out_cpd_als.Algorithm.fval;
    iter_cpd_als=out_cpd_als.Algorithm.iterations;
    
    iter_cpd_als_result(i)=iter_cpd_als;
    cput_cpd_als_result(i)=cput_cpd_als;
    res_cpd_als_result(i)=norm(lab_res(:));
    rel_error_cpd_als_result(i)=norm(lab_res(:))/norm(X(:));

    
    %%%%%%%%%%%%%%%%%%%%%%%%% cpd_nls method
    options_nls = struct;
    options_nls.Compression = false;
    options_nls.Algorithm = @cpd_nls;
    options_nls.AlgorithmOptions.MaxIter = 1500;      
    
    options_nls.AlgorithmOptions.CGMaxIter = 500;      
    timect_rnd = cputime;
    
    [Uhat,out_cpd_nls]=cpd(X_noise,Uinit,options_nls);
    cput_cpd_nls = cputime - timect_rnd;
    lab_res_rand=X-cpdgen(Uhat);
    fval_cpd_nls=out_cpd_nls.Algorithm.fval;
    iter_cpd_nls=out_cpd_nls.Algorithm.iterations;

    iter_cpd_nls_result(i)=iter_cpd_nls;
    cput_cpd_nls_result(i)=cput_cpd_nls;
    res_cpd_nls_result(i)=norm(lab_res_rand(:));
    rel_error_cpd_nls_result(i)=norm(lab_res_rand(:))/norm(X(:));

    Rel_err(count,:)=[rel_error_cp_als_result(i),rel_error_cpd_als_result(i),rel_error_cpd_nls_result(i),rel_error_vp_result(i)];
    Time(count,:)=[cput_als,cput_cpd_als,cput_cpd_nls,cput_vp_0];
    count=count+1;

    if r==15
        fprintf("============================= iter=%d ===================\n",i);
        fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f, res_err is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)),norm(Xres0(:))/norm(X(:)));
        fprintf("iter=%d, cput=%3.4f, the residual of the als method is %3.10f, res_err is %3.10f\n",iter_als,cput_als,norm(als_res(:)),norm(als_res(:))/norm(X(:)));
        fprintf("iter=%d, cput=%3.4f, the residual of the cpd_als method is %3.10f, res_err is %3.10f\n",iter_cpd_als,cput_cpd_als,norm(lab_res(:)), norm(lab_res(:))/norm(X(:)));
        fprintf("iter=%d, cput=%3.4f, the residual of the cpd_nls method is %3.10f, res_err is %3.10f\n",iter_cpd_nls,cput_cpd_nls,norm(lab_res_rand(:)),norm(lab_res_rand(:))/norm(X(:)));
    end
end

fprintf("========================================== rank %d approximation ==========================================\n",r);
fprintf("%d &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f\n",r,mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result));

end 


formatSpec = '%.4f %.4f %.4f %.4f\n';
r_id = fopen('rel_err.txt', 'w');
t_id = fopen('time.txt', 'w');
fprintf(r_id, formatSpec, Rel_err');
fprintf(t_id, formatSpec, Time');

fclose(r_id); 
fclose(t_id); 