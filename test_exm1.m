clc; clear all; close all;


% addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
% addpath("E:\\matlab-code\\variable-projection\\tensorlab")

root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, 'tensor_toolbox')));
addpath(genpath(fullfile(root, 'tensorlab')));
addpath(genpath(fullfile(root, 'utils')));

% case 1
m=100;
n=100;
p=100;

% case 2
% m=30;
% n=40;
% p=1000;

% case 3
% m=100;
% n=100;
% p=10000;


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


    %%%%%%%%%%% test data
    r=12;
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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% vp_pGN method (ours)
    vp_paras.gamma=0;
    vp_paras.itmax=1500;
    vp_paras.xtol=1e-6; 
    vp_paras.ftol=1e-12; 
    
    cg_paras.itmax=500;
    cg_paras.tol=1e-6;
    
    runhist_vp_0 = vp_pGN(U,V,U_true,V_true,W_true,vp_paras,cg_paras);
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
    [P,U0,out] = cp_als(tensor(X),r,...
        'printitn',0 , ...
        'init', Uinit, ...
        'maxiters', 1500, ...
        'tol', 1e-6);

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
    timect = cputime;
    [Xhat,out_cpd_als]=cpd(X,Uinit,options_als);
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
    
    [Uhat,out_cpd_nls]=cpd(X,Uinit,options_nls);
    cput_cpd_nls = cputime - timect_rnd;
    lab_res_rand=X-cpdgen(Uhat);
    fval_cpd_nls=out_cpd_nls.Algorithm.fval;
    iter_cpd_nls=out_cpd_nls.Algorithm.iterations;

    iter_cpd_nls_result(i)=iter_cpd_nls;
    cput_cpd_nls_result(i)=cput_cpd_nls;
    res_cpd_nls_result(i)=norm(lab_res_rand(:));
    rel_error_cpd_nls_result(i)=norm(lab_res_rand(:))/norm(X(:));

    Rel_err(count,:)=[rel_error_cp_als_result(i),rel_error_cpd_als_result(i),rel_error_cpd_nls_result(i),rel_error_vp_result(i)];
    T(count,:)=[cput_als,cput_cpd_als,cput_cpd_nls,cput_vp_0];
    count=count+1;

    fprintf("============================= iter=%d ===================\n",i);
    fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f, res_err is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)),norm(Xres0(:))/norm(X(:)));
    fprintf("iter=%d, cput=%3.4f, the residual of the als method is %3.10f, res_err is %3.10f\n",iter_als,cput_als,norm(als_res(:)),norm(als_res(:))/norm(X(:)));
    fprintf("iter=%d, cput=%3.4f, the residual of the cpd_als method is %3.10f, res_err is %3.10f\n",iter_cpd_als,cput_cpd_als,norm(lab_res(:)), norm(lab_res(:))/norm(X(:)));
    fprintf("iter=%d, cput=%3.4f, the residual of the cpd_nls method is %3.10f, res_err is %3.10f\n",iter_cpd_nls,cput_cpd_nls,norm(lab_res_rand(:)),norm(lab_res_rand(:))/norm(X(:)));
end

fprintf("========================================== rank %d approximation ==========================================\n",r);
fprintf("%d &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f\n",r,mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result));


fprintf("relative error of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(rel_error_cp_als_result),mean(rel_error_cpd_als_result),mean(rel_error_cpd_nls_result),mean(rel_error_vp_result));
fprintf("mean value of residul value of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(res_cp_als_result),mean(res_cpd_als_result),mean(res_cpd_nls_result),mean(res_vp_result));
fprintf("mean value of iteration of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(iter_cp_als_result),mean(iter_cpd_als_result),mean(iter_cpd_nls_result),mean(iter_vp_result));
fprintf("mean value of cpu time of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(cput_cp_als_result),mean(cput_cpd_als_result),mean(cput_cpd_nls_result),mean(cput_vp_result));
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


 