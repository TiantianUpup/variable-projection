clc; clear all; close all;


addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
addpath("E:\\matlab-code\\variable-projection\\tensorlab")


m=50;
n=50;
p=50;

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
Time=zeros(10,4);
count=1;

r_true=3;
for r=2:r_true
for i=1:trial
    % truth factor matrices   
    % colinear factor matrices
    C=col_mat_gen(m,r_true,c);
    U_true=C{1};
    V_true=C{2};
    W_true=C{3};
    T={U_true,V_true,W_true};

    X=generate_cp_tensor(U_true,V_true,W_true);
  

    % %% complete orthogonality
    % U_true = orth(rand(m,r));
    % V_true = orth(rand(n,r));
    % W_true = orth(rand(p,r));

    

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
    relfval_vp = runhist_vp_0.fval;
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
        'tol', 1e-8);

    cput_als = cputime - timect;   
    als_res=full(P)-X;
    relfval_als=out.fval;
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
    [Xhat,out_cpd_als]=cpd(X,Uinit,options_als);
    cput_cpd_als = cputime - timect;
    lab_res=X-cpdgen(Xhat);
    relfval_cpd_als=out_cpd_als.Algorithm.relfval;
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
    relfval_cpd_nls=out_cpd_nls.Algorithm.relfval;
    iter_cpd_nls=out_cpd_nls.Algorithm.iterations;
  
    iter_cpd_nls_result(i)=iter_cpd_nls;
    cput_cpd_nls_result(i)=cput_cpd_nls;
    res_cpd_nls_result(i)=norm(lab_res_rand(:));
    rel_error_cpd_nls_result(i)=norm(lab_res_rand(:))/norm(X(:));

    Rel_err(count,:)=[rel_error_cp_als_result(i),rel_error_cpd_als_result(i),rel_error_cpd_nls_result(i),rel_error_vp_result(i)];
    Time(count,:)=[cput_als,cput_cpd_als,cput_cpd_nls,cput_vp_0];
    count=count+1;

    % fprintf("============================= iter=%d ===================\n",i);
    % fprintf("iter=%d, cput=%3.4f, the residual of X and Xhat for variable_projection_gamma (gamma=0) is %3.20f, res_err is %3.10f\n",iter_vp_0,cput_vp_0,norm(Xres0(:)),norm(Xres0(:))/norm(X(:)));
    % fprintf("iter=%d, cput=%3.4f, the residual of the als method is %3.20f, res_err is %3.10f\n",iter_als,cput_als,norm(als_res(:)),norm(als_res(:))/norm(X(:)));
    % fprintf("iter=%d, cput=%3.4f, the residual of the cpd_als method is %3.20f, res_err is %3.10f\n",iter_cpd_als,cput_cpd_als,norm(lab_res(:)), norm(lab_res(:))/norm(X(:)));
    % fprintf("iter=%d, cput=%3.4f, the residual of the cpd_nls method is %3.20f, res_err is %3.10f\n",iter_cpd_nls,cput_cpd_nls,norm(lab_res_rand(:)),norm(lab_res_rand(:))/norm(X(:)));
end

%fprintf("========================================== rank %d approximation ==========================================\n",r);
fprintf("%d &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f &%3.1f &%.4f &%3.2f\n",r,mean(iter_cp_als_result),mean(rel_error_cp_als_result),mean(cput_cp_als_result),mean(iter_cpd_als_result),mean(rel_error_cpd_als_result),mean(cput_cpd_als_result),mean(iter_cpd_nls_result),mean(rel_error_cpd_nls_result),mean(cput_cpd_nls_result), mean(iter_vp_result),mean(rel_error_vp_result),mean(cput_vp_result));
end

% fprintf("relative error of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(rel_error_cp_als_result),mean(rel_error_cpd_als_result),mean(rel_error_cpd_nls_result),mean(rel_error_vp_result));
% fprintf("mean value of residul value of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f,vp=%.4f\n",mean(res_cp_als_result),mean(res_cpd_als_result),mean(res_cpd_nls_result),mean(res_vp_result));
% fprintf("mean value of iteration of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(iter_cp_als_result),mean(iter_cpd_als_result),mean(iter_cpd_nls_result),mean(iter_vp_result));
% fprintf("mean value of cpu time of cp_als=%.4f, cpd_als=%.4f, cpd_nls=%.4f, vp=%.4f\n",mean(cput_cp_als_result),mean(cput_cpd_als_result),mean(cput_cpd_nls_result),mean(cput_vp_result));
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
% r_id = fopen('rel_err.txt', 'w');
% t_id = fopen('time.txt', 'w');
% fprintf(r_id, formatSpec, Rel_err');
% fprintf(t_id, formatSpec, Time');
% 
% fclose(r_id); 
% fclose(t_id); 
% disp("Time matrix is:");
% disp(Time);

%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
fig=figure();
% semilogy(fval_cpd_als,'-o', 'LineWidth', 1.5); hold all; %'g--s',
% semilogy(fval_als,'-*', 'LineWidth', 1.5); hold all; % 'g--s',
% semilogy(fval_cpd_nls, '-square','LineWidth', 1.5); hold all; % 'b-*',
% semilogy(fval_vp, '-x','LineWidth', 1.5); 


% plot(log10(fval_als),'-*', 'LineWidth', 1.5); hold all; % 'g--s',
% plot(log10(fval_cpd_als),'-o', 'LineWidth', 1.5); hold all; %'g--s',
% plot(log10(fval_cpd_nls), '-square','LineWidth', 1.5); hold all; % 'b-*',
% plot(log10(fval_vp), '-x','LineWidth', 1.5); 



plot(log10(relfval_als),'-*', 'LineWidth', 1.5); hold all; % 'g--s',
plot(log10(relfval_cpd_als),'-o', 'LineWidth', 1.5); hold all; %'g--s',
plot(log10(relfval_cpd_nls), '-square','LineWidth', 1.5); hold all; % 'b-*',
plot(log10(relfval_vp), '-x','LineWidth', 1.5); 

% semilogy(fval_vp); hold all;
% semilogy(fval_cpd_als); hold all;
% semilogy(fval_cpd_nls);
ylabel('log_{10}(relfval)','FontSize',20,'FontWeight','bold'); 
xlabel('iterations','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',13,'FontWeight','bold');
%title('Convergence plot'); 


% max_iter = max([length(fval_cpd_als), length(fval_als), length(fval_cpd_nls), length(fval_vp)]);
% xlim([-10, max_iter]);

% x_range = xlim;
% xlim([-15, x_range(2)]);
% 
% y_range = ylim;
% ylim([-5, y_range(2)]);
%%xticks(0:50:max_iter);
legend('cp\_als','cpd\_als','cpd\_nls','vp\_pGN','FontSize', 15,'FontWeight', 'bold','Location', 'southeast');  % 
xlim([0,50]);
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

print(fig,'col-3','-dpdf') 
 