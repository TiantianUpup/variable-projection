clc; clear all; close all;

addpath('../compared_method/tensorlab');
addpath('../compared_method/tensor_toolbox');
addpath('../utils');
addpath('../algo');




% m=30;
% n=40;
% p=1000;

% m=15;
% n=15;
% p=15;

m=50;
n=50;
p=50;

% % m=30;
% % n=30;
% % p=40;

% m=50;
% n=50;
% p=50;

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


% U_true=rand(m,r);
% V_true=rand(n,r);
% W_true=rand(p,r);
Rel_err=zeros(10,4);
Time=zeros(10,4);
count=1;

l1=5;
l2=5;
for r=1:16
for i=1:trial
    % truth factor matrices
    r_true=15;
    U_true=rand(m,r_true);
    V_true=rand(n,r_true);
    W_true=rand(p,r_true);

    %MA=khatrirao(U_true,V_true)*W_true';
    X=generate_cp_tensor(U_true,V_true,W_true);

    % N=randn(m,n,p);
    % rho_t=0.0005*(norm(X(:))/norm(N(:)));
    % X_noise=X+rho_t*N;

    % mn=m*n;
    % N=randn(mn,p);
    % rho=0.1*((norm(MA(:))/norm(N(:))));
    % MA_noise=MA+rho*N;

    %U_n=randn(m,r_true);
    %rho_u=0.005*(norm(U_true(:))/norm(U_n(:)));
    %rho_u=(100/l1-1)^(-1/2)*(norm(U_true(:))/norm(U_n(:)));
    %disp(rho_u);
    %U_noise=U_true+rho_u*U_n;

    % W_n=randn(p,r_true);
    % rho_w=0.01*(norm(W_true(:))/norm(W_n(:)));
    % % rho_u=(100/l1-1)^(-1/2)*(norm(U_true(:))/norm(U_n(:)));
    % % disp(rho_u);
    % W_noise=W_true+rho_w*W_n;

    % rho = 0.02;                 
    % alpha = 10*max(abs(W_true(:)));  
    % W_noise = W_true;
    % mask = rand(size(W_true)) < rho;
    % W_noise(mask) = W_true(mask) + alpha*randn(nnz(mask),1);

    % tau=0.01*norm(W_true,'fro')/sqrt(numel(W_true));
    % E=randn(size(W_true));
    % W_noise=W_true+tau*E;

    % rho = 0.05;
    % W_noise = W_true;
    % mask = rand(size(W_true)) < rho;
    % sgn = sign(rand(nnz(mask),1)-0.5);
    % W_noise(mask) = 100*sgn;

    [U,S,V] = svd(W_true,'econ');
    E = -S(end,end)*U(:,end)*V(:,end)';
    W_true = W_true + E;
    %fprintf("===============s the rank of W_noise is %d =============\n", rank(W_noise));

    T={U_true,V_true,W_true};
    X_noise=generate_cp_tensor(U_true,V_true,W_true);

    % [U_noise,N]=noisy(U_true,10);
    % T={U_noise,V_true,W_true};
    % X_noise=generate_cp_tensor(U_noise,V_true,W_true);

    % [W_noise,N]=noisy(W_true,5);
    % T={U_true,V_true,W_noise};
    % X_noise=generate_cp_tensor(U_true,V_true,W_noise);

    % U_n_1=randn(m,r_true);
    % U_n_temp=U_n_1.*U_noise;
    % U_noise_f=U_noise+(100/l2-1)^(-1/2)*(norm(U_noise(:))/norm(U_n_temp(:)))*U_n_temp;

    
    %V_n=randn(n,r_true);
    % rho_v=0.005*(norm(V_true(:))/norm(V_n(:)));
    % %rho_v=(100/l1-1)^(-1/2)*(norm(V_true(:))/norm(V_n(:)));
    % %disp(rho_v);
    % V_noise=V_true+rho_v*V_n;

    % V_n_1=randn(n,r_true);
    % V_n_temp=V_n_1.*V_noise;
    % V_noise_f=V_noise+(100/l2-1)^(-1/2)*(norm(V_noise(:))/norm(V_n_temp(:)))*V_n_temp;

    % T={U_noise,V_noise,W_true};
    % X_noise=generate_cp_tensor(U_noise,V_noise,W_true);

    % T={U_noise_f,V_noise_f,W_true};
    % X_noise=generate_cp_tensor(U_noise_f,V_noise_f,W_true);


    % U_noise=U_true;
    % V_noise=V_true;
    % T=noisy({U_true,V_true,W_true},60);

    %X_noise=generate_cp_tensor(U_noise,V_noise,W_noise);
    % U_noise=T{1};
    % V_noise=T{2};
    % W_noise=T{3};
    % X_noise=generate_cp_tensor(U_noise,V_noise,W_noise);

    % [W_noise,N]=noisy(W_true,25);

    % sigma_t=20*log10(norm(W_true(:))/norm(N(:)));
    % fprintf("======================== sigma is %3.10f ===================\n",sigma_t);

    % W_n=randn(p,r_true);
    % rho_w=0.005*(norm(W_true(:))/norm(W_n(:)));
    % % disp(rho_w);
    % W_noise=W_true+rho_w*W_n;
    % T={U_true,V_true,W_noise};
    % X_noise=generate_cp_tensor(U_true,V_true,W_noise);

    %%%%%%%%%%% test data
    % r=8;
    paras.m=m;
    paras.n=n;   
    paras.r=r;

    % U = orth(rand(m,r));
    % V = orth(rand(n,r));

    % initial point (random)
    U = rand(m,r);
    U = proj_oblique(U);
    V = rand(n,r);
    V = proj_oblique(V);
    W = rand(p,r);
    W = proj_oblique(W);
    Uinit={U,V,W};
    % Uinit_vp={U,V};

    % gevd initialization
    % cput_gevd=cputime;
    % [Uinit,output] = cpd_gevd(X,r);
    % fprintf("--------------------------- cpd_gevd costs %3.10f ----------------------------\n",cputime-cput_gevd);

    % U=Uinit{1};
    % V=Uinit{2};
    
    % fprintf("the size of U is\n");
    % disp(size(U));

    vp_paras.itmax=1500;
    vp_paras.xtol=1e-6; 
    vp_paras.ftol=1e-12; 
    
    cg_paras.itmax=500;
    cg_paras.tol=1e-6;

    %profile on;
    runhist_vp_0 = vp_pGN(T, Uinit, vp_paras, cg_paras);
    
    %profile off;
    %profile viewer;
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

    Time(count,:)=[cput_als,cput_cpd_als,cput_cpd_nls,cput_vp_0];
    Rel_err(count,:)=[rel_error_cp_als_result(i),rel_error_cpd_als_result(i),rel_error_cpd_nls_result(i),rel_error_vp_result(i)];
    count=count+1;

    % fprintf("============================= iter=%d ===================\n",i);
    % fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.10f, res_err is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)),norm(Xres0(:))/norm(X(:)));
    % fprintf("iter=%d, cput=%3.4f, the residual of the als method is %3.10f, res_err is %3.10f\n",iter_als,cput_als,norm(als_res(:)),norm(als_res(:))/norm(X(:)));
    % fprintf("iter=%d, cput=%3.4f, the residual of the cpd_als method is %3.10f, res_err is %3.10f\n",iter_cpd_als,cput_cpd_als,norm(lab_res(:)), norm(lab_res(:))/norm(X(:)));
    % fprintf("iter=%d, cput=%3.4f, the residual of the cpd_nls method is %3.10f, res_err is %3.10f\n",iter_cpd_nls,cput_cpd_nls,norm(lab_res_rand(:)),norm(lab_res_rand(:))/norm(X(:)));
end

fprintf("========================================== rank %d approximation ==========================================\n",r);
fprintf("%d &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f\n",r,mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result));

% fprintf("relative error of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(rel_error_cp_als_result),mean(rel_error_cpd_als_result),mean(rel_error_cpd_nls_result),mean(rel_error_vp_result));
% fprintf("mean value of residul value of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(res_cp_als_result),mean(res_cpd_als_result),mean(res_cpd_nls_result),mean(res_vp_result));
% fprintf("mean value of iteration of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(iter_cp_als_result),mean(iter_cpd_als_result),mean(iter_cpd_nls_result),mean(iter_vp_result));
% fprintf("mean value of cpu time of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(cput_cp_als_result),mean(cput_cpd_als_result),mean(cput_cpd_nls_result),mean(cput_vp_result));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf("vp, min=%3.4f, max=%3.4f\n",min(rel_error_vp_result),max(rel_error_vp_result));
% fprintf("cp_als, min=%3.4f, max=%3.4f\n",min(rel_error_cp_als_result),max(rel_error_cp_als_result));
% fprintf("cpd_als, min=%3.4f, max=%3.4f\n",min(rel_error_cpd_als_result),max(rel_error_cpd_als_result));
% fprintf("cpd_nls, min=%3.4f, max=%3.4f\n",min(rel_error_cpd_nls_result),max(rel_error_cpd_nls_result));

% Result = [res_cp_als_result-res_vp_result, res_cpd_als_result-res_vp_result, res_cpd_nls_result-res_vp_result];
% prob=sum(all(Result >= 0, 2));
% % prob=
% max_gap=max(max(Result(Result >= 0)));
% fprintf("the probability is %3.1f, the max residual is %3.2f\n",prob/trial,max_gap);
% disp("Result matrix is:");
% disp(Result);

% disp("Rel_err matrix is:");
% disp(Rel_err);

% formatSpec = '%.4f %.4f %.4f %.4f\n';
% formatSpec_t = '%.4f %.4f %.4f %.4f\n';
% r_id = fopen('rel_err.txt', 'w');
% t_id = fopen('time.txt', 'w');
% fprintf(r_id, formatSpec, Rel_err');
% fprintf(t_id, formatSpec_t, Time');

% fclose(r_id); 
% fclose(t_id); 
% disp("Time matrix is:");
% disp(Time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% fprintf(r, mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result),mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result));

% Result_res = [res_cp_als_result-res_vp_result, res_cpd_als_result-res_vp_result, res_cpd_nls_result-res_vp_result];
% prob_res=sum(all(Result_res >= 0, 2));
% % prob=
% max_gap=max(max(Result_res(Result_res >= 0)));
% fprintf("the probability is %3.1f, the max residual is %3.2f\n",prob_res/trial,max_gap);
% disp("Result_res matrix is:");
% disp(Result_res);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
if ismember(r, [5, 10, 15])
    fig=figure();
    % semilogy(fval_cpd_als,'-o', 'LineWidth', 1.5); hold all; %'g--s',
    % semilogy(fval_als,'-*', 'LineWidth', 1.5); hold all; % 'g--s',
    % semilogy(fval_cpd_nls, '-square','LineWidth', 1.5); hold all; % 'b-*',
    % semilogy(fval_vp, '-x','LineWidth', 1.5); 

    plot(log10(fval_cpd_als),'-o', 'LineWidth', 1.5); hold all; %'g--s',
    plot(log10(fval_als),'-*', 'LineWidth', 1.5); hold all; % 'g--s',
    plot(log10(fval_cpd_nls), '-square','LineWidth', 1.5); hold all; % 'b-*',
    plot(log10(fval_vp), '-x','LineWidth', 1.5); 

    % semilogy(fval_vp); hold all;
    % semilogy(fval_cpd_als); hold all;
    % semilogy(fval_cpd_nls);
    ylabel('objective function','FontSize',20,'FontWeight','bold'); 
    xlabel('iterations','FontSize',20,'FontWeight','bold');
    set(gca,'FontSize',13,'FontWeight','bold');
    %title('Convergence plot'); 


    % max_iter = max([length(fval_cpd_als), length(fval_als), length(fval_cpd_nls), length(fval_vp)]);
    % xlim([-10, max_iter]);

    x_range = xlim;
    xlim([-15, x_range(2)]);

    % y_range = ylim;
    % ylim([-5, y_range(2)]);
    %%xticks(0:50:max_iter);
    legend('cpd\_als','cp\_als','cpd\_nls','vp\_pGN','FontSize', 15,'FontWeight', 'bold');  % 
    hold off;

    %grid on;  

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

    filename = sprintf('noise-%d', r); 
    print(fig, filename, '-dpdf');
end