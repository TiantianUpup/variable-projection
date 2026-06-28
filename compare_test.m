clc
clear all
close all


addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
addpath("E:\\matlab-code\\variable-projection\\tensorlab")
% m=10;
% n=10;
% p=100;

m=30;
n=40;
p=40;

% m=2; 
% n=3;
% p=3;

r=15;
% % truth factor matrices
U_true=rand(m,r);


%and(m,1);
% for i=1:r
%     U_true(:,i)=u+0.0001*randn(m,1);
% end 


%V_true = 2*U_true;
V_true=rand(n,r);
% %V=2*U;
% V_true=zeros(n,r);
% v=rand(n,1);
% for i=1:r
%     V_true(:,i)=v+0.0001*randn(n,1);
% end 

%W_true = 3*U_true;
W_true=rand(p,r);

% %W=3*U;
% 
%ill-conditioned matrix 
% kappa=1e15;
% Z = orth(rand(p, r));  
% R = orth(rand(r, r)); 
    
% s = logspace(0, -log10(kappa), r);

% W_true = Z*diag(s)*R';

% W_true=zeros(p,r);
% w=rand(p,1);
% kappa=1e-6;
% for i=1:r
%     W_true(:,i)=kappa*rand(p,1);
%     kappa=kappa*10;
% end    

% fprintf("------------------------------ the condition of the W is %3.10f ------------------\n",cond(W_true));


% % truth factor matrices
% U = randn(m,r);
% V = randn(n,r);
% W = randn(p,r);

%% complete orthogonality
% U = orth(rand(m,r));
% V = orth(rand(n,r));
% W = orth(rand(p,r));

X=generate_cp_tensor(U_true,V_true,W_true);
% rho=0.00005; case 1
%rho=0.005;
N=randn(m,n,p);
rho_t=0.1*(norm(X(:))/norm(N(:)));
X_noise=X+rho_t*N;

T_size=sum(size(X_noise));
fprintf("---------------- size(X_noise) is %d ----------------------\n",T_size);

% A concrete example
% X=zeros(m,n,p);
% U_true=rand(m,2);
% U_1=U_true(:,1);
% U_2=U_true(:,2);
% V_true=rand(n,2);
% V_1=V_true(:,1);
% V_2=V_true(:,2);
% W_true=rand(p,2);
% W_1=W_true(:,1);
% W_2=W_true(:,2);
% X=reshape(kron(kron(W_2,V_1),U_1),[m,n,p])+reshape(kron(kron(W_1,V_2),U_1),[m,n,p])+reshape(kron(kron(W_1,V_1),U_2),[m,n,p]);

% U_true=rand(m,2);
% U_1=U_true(:,1);
% U_2=U_true(:,2);
% V_true=rand(n,2);
% V_1=V_true(:,1);
% V_2=V_true(:,2);
% W_true=rand(p,2);
% W_1=W_true(:,1);
% W_2=W_true(:,2);
% X=reshape(kron(kron(W_1,V_2),U_1),[m,n,p])+reshape(kron(kron(W_1,V_1),U_2),[m,n,p])+reshape(kron(kron(W_2,V_1),U_1),[m,n,p]);

%%%%%%%%%%% test data
r=15;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                ; % rank-r approximation

MA=khatrirao(U_true,V_true)*W_true';

% noisy case
mn=m*n;
N=randn(mn,p);
rho=0.1*((norm(MA(:))/norm(N(:))));
MA_noise=MA+rho*N;


paras.m=m;
paras.n=n;
paras.r=r;

% U = orth(rand(m,r));
% V = orth(rand(n,r));

U = rand(m,r);
V = rand(n,r);

 
%%%%%%%%%%% variable projection method for original problem
% runhist_ori = gauss_newton_original(U,V,MA, paras);
% Uhat = runhist_ori.U;
% Vhat = runhist_ori.V;
% cput_ori = runhist_ori.cput;
% iter_ori = runhist_ori.iter;
% 
% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% 
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% Xres=Xhat-X;

%%%%%%%%%%% variable projection method for gamma=0

vp_paras.gamma=0;
vp_paras.itmax=1500;
vp_paras.tol=1e-6; 
vp_paras.mu=0;
 
cg_paras.itmax=500;
cg_paras.tol=1e-6;
cg_paras.lambda=0;

%profile on;
%runhist_vp_0 = variable_projection_simple(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
runhist_vp_0 = variable_projection_simple(U,V,U_true,V_true,W_true,MA_noise,vp_paras,cg_paras,paras);
%profile off;
%profile viewer;
%runhist_vp_0 = variable_projection_gamma_acg(U,V,MA,vp_paras,cg_paras,paras);
%runhist_vp_0 = variable_projection_simple(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
Uhat_1 = runhist_vp_0.U;
Vhat_1 = runhist_vp_0.V;
What_1 = runhist_vp_0.W;
fval_vp=runhist_vp_0.fval;
cput_vp_0 = runhist_vp_0.cput;
iter_vp_0 = runhist_vp_0.iter;

Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
Xres0=Xhat_1-X;


% runhist_vp_0 = variable_projection_fopt(U,V,MA,vp_paras,cg_paras,paras);
% Uhat_1 = runhist_vp_0.U;
% Vhat_1 = runhist_vp_0.V;
% What_1 = runhist_vp_0.W;
% cput_fopt = runhist_vp_0.cput;
% iter_fopt = runhist_vp_0.iter;
% 
% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_fopt=Xhat_1-X;


% runhist_mcg = variable_projection_simple_acg(U,V,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras);
runhist_mcg = variable_projection_simple_acg(U,V,U_true,V_true,W_true,MA_noise,vp_paras,cg_paras,paras);
Uhat_1 = runhist_mcg.U;
Vhat_1 = runhist_mcg.V;
What_1 = runhist_mcg.W;
cput_mcg = runhist_mcg.cput;
iter_mcg = runhist_mcg.iter;

Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
Xres_mcg=Xhat_1-X;

% profile on;
% runhist_pcg = variable_projection_pcg(U,V,MA,vp_paras,cg_paras,paras);
% Uhat_1 = runhist_pcg.U;
% Vhat_1 = runhist_pcg.V;
% What_1 = runhist_pcg.W;
% cput_pcg = runhist_pcg.cput;
% iter_pcg = runhist_pcg.iter;

% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_pcg=Xhat_1-X;

% profile off;
% profile viewer;

% runhist_vp_opt = variable_projection_opt(U,V,MA,vp_paras,cg_paras,paras);
% Uhat_1 = runhist_vp_opt.U;
% Vhat_1 = runhist_vp_opt.V;
% cput_vp_opt = runhist_vp_opt.cput;
% iter_vp_opt = runhist_vp_opt.iter;
% 
% kr=khatrirao(Uhat_1,Vhat_1);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What_1=MA'*krmp';
% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_opt=Xhat_1-X;


% 
% % BFGS method
% runhist_bfgs = BFGS(U,V,MA,paras);
% Uhat_1 = runhist_bfgs.U;
% Vhat_1 = runhist_bfgs.V;
% cput_bfgs = runhist_bfgs.cput;
% iter_bfgs = runhist_bfgs.iter;
% 
% kr=khatrirao(Uhat_1,Vhat_1);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What_1=MA'*krmp';
% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_bfgs=Xhat_1-X;



% %%%%%%%%%%% variable projection method for gamma=1
% vp_paras.gamma=1;

% runhist_vp_1 = variable_projection_gamma(U,V,MA,vp_paras,cg_paras,paras);
% Uhat = runhist_vp_1.U;
% Vhat = runhist_vp_1.V;
% cput_vp_1 = runhist_vp_1.cput;
% iter_vp_1 = runhist_vp_1.iter;

% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat_gamma=generate_cp_tensor(Uhat,Vhat,What);
% Xres1=Xhat_gamma-X;


% % %%%%%%%%%%% gradient descent method
% alpha=1e-4;
% beta=0.8;
% epsilon=1e-6;
% runhist_gd=gradient_descent(U,V,MA,alpha,beta,epsilon);
% %runhist_gd=agd(U,V,MA,alpha,beta,epsilon);
% Uhat = runhist_gd.U;
% Vhat = runhist_gd.V;
% cput = runhist_gd.cput;
% iter = runhist_gd.iter;
% 
% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% resg=X-Xhat;

%%%%%%%%%%% variable projection method for gamma=1
% gamma=1;
% runhist_vp = variable_projection_gamma(U,V,MA,gamma,paras);
% Uhat = runhist_vp.U;
% Vhat = runhist_vp.V;
% cput_vp_1 = runhist_vp.cput;
% iter_vp_1 = runhist_vp.iter;
% 
% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% Xreso=Xhat-X;

% %%%%%%%%%%% two-phase method
% epsilon=1e-2;
% runhist_gd_gn = gd_gn_method(U,V,MA,alpha,beta,epsilon);
% Uhat = runhist_gd_gn.U;
% Vhat = runhist_gd_gn.V;

% kr=khatrirao(Uhat,Vhat);
% krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
% What=MA'*krmp';
% Xhat=generate_cp_tensor(Uhat,Vhat,What);
% Xres_gd_gn=Xhat-X;
% iter_gd=runhist_gd_gn.iter_gd;
% iter_gn=runhist_gd_gn.iter_gn;
% cput_gd_gn=runhist_gd_gn.cput;
% cput_gn=runhist_gd_gn.cput_gn;
% cput_gd=runhist_gd_gn.cput_gd;

% runhist_two_phase = two_phase_als(X,MA,vp_paras,cg_paras,paras);
% Uhat_1 = runhist_two_phase.U;
% Vhat_1 = runhist_two_phase.V;
% What_1 = runhist_two_phase.W;
% cput_vp = runhist_two_phase.cput_vp;
% iter_vp = runhist_two_phase.iter_vp;
% cput_als_2 = runhist_two_phase.cput_als;
% iter_als_2 = runhist_two_phase.iter_als;

% Xhat_1=generate_cp_tensor(Uhat_1,Vhat_1,What_1);
% Xres_2=Xhat_1-X;

% tensor_toolbox method
timect = cputime;
[P,U0,out] = cp_als(tensor(X_noise),r,...
    'maxiters', 1500, ...
    'tol', 1e-6);
% [P,U0,out] = cp_als(tensor(X),r,...
%     'maxiters', 1500, ...
%     'tol', 1e-6);
cput_als = cputime - timect;   
als_res=full(P)-X;
fval_als=out.fval;
iter_als=out.iters;

% tensorlab method
options_als = struct;
%options_als.Compression = true;
options_als.Algorithm = @cpd_als;
options_als.AlgorithmOptions.MaxIter = 1500;      % Default 500
options_als.AlgorithmOptions.CGMaxIter = 500;     % Default 15
options_als.Initialization = @cpd_rnd;
timect = cputime;
[Xhat,out_cpd_als]=cpd(X_noise,r,options_als);
%[Xhat,out_cpd_als]=cpd(X,r,options_als);
cput_tensorlab = cputime - timect;
lab_res=X-cpdgen(Xhat);
fval_cpd_als=out_cpd_als.Algorithm.fval;
iter_cpd_als=out_cpd_als.Algorithm.iterations;
disp(out_cpd_als.Algorithm);

% nls method
options_nls = struct;
%options_nls.Compression = true;
options_nls.Algorithm = @cpd_nls;
options_nls.AlgorithmOptions.MaxIter = 1500;      % Default 500
options_nls.AlgorithmOptions.CGMaxIter = 15;     % Default 15
% timect = cputime;
% [Xhat,out]=cpd(X,r,options_nls);
% cput_nls = cputime - timect;
% lab_nls=X-cpdgen(Xhat);
% iter_cpd_nls=out.Algorithm.iterations;
% %disp(out);
% disp(out.Algorithm);


% nls method with rand initial point
options_nls.Initialization = @cpd_rnd;
options_nls.Compression = false;
timect_rnd = cputime;
[Uhat,out_cpd_nls]=cpd(X_noise,r,options_nls);
%[Uhat,out_cpd_nls]=cpd(X,r,options_nls);
cput_rand = cputime - timect_rnd;
lab_res_rand=X-cpdgen(Uhat);
fval_cpd_nls=out_cpd_nls.Algorithm.fval;
iter_cpd_rand=out_cpd_nls.Algorithm.iterations;
disp(out_cpd_nls.Algorithm);
disp(out_cpd_nls.Compression);
disp(out_cpd_nls.Preprocessing);
% disp("------------------- delta is -----------------------");
% disp(out_cpd_nls.Algorithm.delta);
disp("===============================================");


% timect_rnd = cputime;
% options_nls.Compression = false;
% [Uhat,out_cpd_nls]=cpd(X_noise,r,options_nls);
% %[Uhat,out_cpd_nls]=cpd(X,r,options_nls);
% cput_rand_compress = cputime - timect_rnd;
% fval_cpd_nls_compress=out_cpd_nls.Algorithm.fval;
% iter_cpd_rand_compress=out_cpd_nls.Algorithm.iterations;
% lab_res_rand_compress=X-cpdgen(Uhat);

fprintf("============================= final result ===================\n");
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma is %3.10f\n",iter_vp_1,cput_vp_1, norm(Xreso(:))^2);
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for gauss_newton_original is %3.10f\n",iter_ori,cput_ori,norm(Xres(:))^2);

fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f, res_err is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)),norm(Xres0(:))/norm(X(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_foptt (gamma=0) is %3.10f\n",iter_fopt,cput_fopt,norm(Xres0(:)));
%fprintf("iter_vp=%d, iter_als=%d, cput_vp=%3.4f, cput_als=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f\n",iter_vp,iter_als_2,cput_vp,cput_als_2,norm(Xres_2(:)));

fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_mcg (gamma=0) is %3.10f, res_err is %3.10f\n",iter_mcg,cput_mcg,norm(Xres_mcg(:)),norm(Xres_mcg(:))/norm(X(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_pcg (gamma=0) is %3.10f\n",iter_pcg,cput_pcg,norm(Xres_pcg(:)));


%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_opt is %3.10f\n",iter_vp_opt,cput_vp_opt,norm(Xres_opt(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for bfgs is %3.10f\n",iter_bfgs,cput_bfgs,norm(Xres_bfgs(:)));

%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=1) is %3.10f\n",iter_vp_1,cput_vp_1,norm(Xres1(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for gradient_descent is %3.10f\n",iter,cput,norm(resg(:)));
%fprintf("gd_iter=%d,gn_iter=%d,cput_gn=%3.4f,cput_gd=%3.4f,cput=%3.4f, the residual of X and Xhat for gradient_descent is %3.10f\n",iter_gd,iter_gn,cput_gn,cput_gd,cput_gd_gn,norm(Xres_gd_gn(:)));
%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm,cput_lm,norm(Xres_lm(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_1,cput_lm_1,norm(Xres_lm_1(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_2,cput_lm_2,norm(Xres_lm_2(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_3,cput_lm_3,norm(Xres_lm_3(:)));
% fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for LM method is %3.10f\n",iter_lm_4,cput_lm_4,norm(Xres_lm_4(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu1(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu2(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu3(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu4(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu5(:)));
% fprintf("the residual of X and Xhat for LM method is %3.10f\n",norm(Xres_mu6(:)));

%fprintf("iter=%d,cput=%3.4f, the residual of X and Xhat for hybrid method is %3.10f\n",iter_hy,cput_hy,norm(Xres_hy(:)));

fprintf("iter=%d,cput=%3.4f,the residual of the als method is %3.10f, res_err is %3.10f\n",iter_als,cput_als,norm(als_res(:)),norm(als_res(:))/norm(X(:)));
fprintf("iter=%d, cput=%3.4f,the residual of the cpd_als method is %3.10f, res_err is %3.10f\n",iter_cpd_als,cput_tensorlab,norm(lab_res(:)), norm(lab_res(:))/norm(X(:)));
% fprintf("iter=%d, cput=%3.4f,the residual of the cpd_nls method is %3.10f\n",iter_cpd_nls,cput_nls,norm(lab_nls(:)));
fprintf("iter=%d, cput=%3.4f,the residual of the cpd_nls method (random initialization) is %3.10f, res_err is %3.10f\n",iter_cpd_rand,cput_rand,norm(lab_res_rand(:)),norm(lab_res_rand(:))/norm(X(:)));
fprintf("iter=%d, cput=%3.4f,the residual of the cpd_nls method (random initialization) is %3.10f, res_err is %3.10f\n",iter_cpd_rand_compress,cput_rand_compress,norm(lab_res_rand_compress(:)),norm(lab_res_rand_compress(:))/norm(X(:)));

fprintf("rho is %3.8f,rho_t is %3.8f\n",rho,rho_t);
%fprintf("&%d &%3.4f &%3.2f &%d &%3.4f &%3.2f &%d &%3.4f &%3.2f &%d &%3.4f &%3.2f\n",iter_mcg,norm(Xres_mcg(:))/norm(X(:)),cput_mcg,iter_als,norm(als_res(:))/norm(X(:)),cput_als,iter_cpd_als,norm(lab_res(:))/norm(X(:)),cput_tensorlab,iter_cpd_rand,norm(lab_res_rand(:))/norm(X(:)),cput_rand);
fprintf("&%d &%3.4f &%3.2f &%d &%3.4f &%3.2f &%d &%3.4f &%3.2f &%d &%3.4f &%3.2f\n",iter_vp_0,norm(Xres0(:))/norm(X(:)),cput_vp_0,iter_als,norm(als_res(:))/norm(X(:)),cput_als,iter_cpd_als,norm(lab_res(:))/norm(X(:)),cput_tensorlab,iter_cpd_rand,norm(lab_res_rand(:))/norm(X(:)),cput_rand);

x1=1:length(fval_vp);
x2=1:length(fval_als);
%fval_cpd_als=out_cpd_als.Algorithm.fval;
x3=1:length(fval_cpd_als);
%fval_cpd_nls=out_cpd_nls.Algorithm.fval;
x4=1:length(fval_cpd_nls);

% figure;
% plot(x1, fval_vp, 'r-o', 'LineWidth', 1.5); hold on;
% plot(x2, fval_cpd_als, 'g--s', 'LineWidth', 1.5);
% hold on;
% plot(x3, fval_cpd_nls, 'b-*', 'LineWidth', 1.5);
% ylabel('Objective function'); xlabel('Iteration');
% title('Convergence plot'); legend('ours','cpd\_als','cpd\_nls')
% hold off;

fig=figure();
semilogy(fval_vp, 'LineWidth', 2); hold all; %'r-o',
semilogy(fval_als, 'LineWidth', 2); hold all; %'g--s',
semilogy(fval_cpd_als, 'LineWidth', 2); hold all; % 'g--s',
semilogy(fval_cpd_nls, 'LineWidth', 2);  % 'b-*',
% semilogy(fval_vp); hold all;
% semilogy(fval_cpd_als); hold all;
% semilogy(fval_cpd_nls);
ylabel('objective function','FontSize',20,'FontWeight','bold'); 
xlabel('iterations','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',13,'FontWeight','bold');
%title('Convergence plot'); 
legend('ours','cp\_als','cpd\_als','cpd\_nls','FontSize', 15,'FontWeight', 'bold');  % 
hold off;

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - 1.1*ti(2) -1.1* ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,'noise-10','-dpdf')

% fprintf("====== U ======\n");
% disp(U);
                                                                                                                                                                                  
% fprintf("====== V ======\n");
% disp(V);

% fprintf("====== W ======\n");
% disp(W);
