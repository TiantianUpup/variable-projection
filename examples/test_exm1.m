%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo_exm1 reproduces the numerical experiments in Subsection 6.1
% (Tables 1-3).
%
% This example considers tensors with a large third dimension p.
%
% The experiments reported in the paper consider
% p = 100, 1000, and 10000 with different values of m and n.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;
addpath('../compared_method/tensor_toolbox');
addpath('../compared_method/tensorlab');
addpath('../utils');
addpath('../algo');



% Case 1: corresponds to Table 1
m=100;
n=100;
p=100;

% Case 2: corresponds to Table 2
% m=30;
% n=40;
% p=1000;

% Case 3: corresponds to Table 3
% m=100;
% n=100;
% p=10000;


rng(24);
trial=10;
r_true=15;

filename = sprintf('../results/ex1-%d.txt', p);
fid = fopen(filename,'w');

myfprintf(fid, "================================ The results obtained by the four methods.=========================\n");
myfprintf(fid, "      ||        cp_als       ||       cpd_als       ||       cpd_nls       ||       vp_pGN        || \n");
myfprintf(fid, " rank ||  rel_er    time     ||  rel_er     time    ||  rel_er     time    ||  rel_er    time     || \n");


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main loop
for r=[12,12,13,14]
for i=1:trial
    % truth factor matrices
    U_true=rand(m,r_true);
    V_true=rand(n,r_true);
    W_true=rand(p,r_true);
    T={U_true,V_true,W_true};

    X=generate_cp_tensor(U_true,V_true,W_true);

    % initial point
    U = rand(m,r);
    U = proj_oblique(U);
    V = rand(n,r);
    V = proj_oblique(V);
    W = rand(p,r);
    W = proj_oblique(W);
    
    Uinit={U,V,W};
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%% vp_pGN method (ours)
    vp_paras.itmax=1500;
    vp_paras.xtol=1e-6; 
    vp_paras.ftol=1e-12; 
    
    cg_paras.itmax=500;
    cg_paras.tol=1e-6;

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
    Time(count,:)=[cput_als,cput_cpd_als,cput_cpd_nls,cput_vp_0];
    count=count+1;
end

myfprintf(fid," %d    ||  %.4f     %3.2f    ||  %.4f      %3.2f   ||  %.4f      %3.2f   ||  %.4f      %3.2f   || \n",r,mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result),mean(rel_error_vp_result),mean(cput_vp_result));

end


r_fn = sprintf('../results/rel_err_exm1-%d.txt', p);
t_fn = sprintf('../results/time_exm1-%d.txt', p);
r_id = fopen(r_fn, 'w');
t_id = fopen(t_fn, 'w');
formatSpec = '%.4f %.4f %.4f %.4f\n';
fprintf(r_id, formatSpec, Rel_err');
fprintf(t_id, formatSpec, Time');

fclose(r_id); 
fclose(t_id); 
fclose(fid); 
 