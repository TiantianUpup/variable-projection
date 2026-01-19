m=100;
n=100;
p=100;
r=10;

addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
U=rand(m,r);
V=rand(n,r);
W=rand(p,r);

UVkr=khatrirao(U,V);
MA=UVkr*W';

% method 1
% tic
% cput_ori=cputime;
% UVkr=khatrirao(U,V);
% Pro=proj_supp(UVkr);
% PMA=Pro*MA;
% %f=fval_opt(PMA);
% f=0.5*norm(PMA(:))^2;
% fprintf("-------------- the ori method costs %3.50f\n", cputime-cput_ori);
% toc


% method 2
tic
cput_opt=cputime;
[P, Sigma, S] = qr(U);
[Q, Tau, T] = qr(V);
S=S';
T=T';

P_tru=P(:,1:r);
%disp(size(P_tru));
Sigma_tru=Sigma(1:r,:);
% disp(size(Sigma_tru));
% disp(size(S));

Q_tru=Q(:,1:r);
Tau_tru=Tau(1:r,:);
X_hat=Sigma_tru*S;
% fprintf("-------------- the size of the X_hat is --------------------------\n");
% disp(size(X_hat));
Y_hat=Tau_tru*T;
XYkr=khatrirao(X_hat,Y_hat);

% fprintf("the size of XYkr is \n");
% disp(size(XYkr));
[P_hat,Sigma_hat,Q_hat]=svd(XYkr,"econ");

% fprintf("the size of P_hat is \n");
% disp(size(P_hat));
PQkr=kron(P_tru,Q_tru);

PQP=PQkr*P_hat;
Proj=PQP*PQP';

R=MA-Proj*MA;
f_opt=0.5*norm(R(:))^2;



cput=cputime;
toc

% fval_res=f-f_opt;

% fprintf("------------------------------ the residual of the fval is %3.40f -----------------------\n",fval_res);

% tic
% PQP=kron(P_tru,Q_tru)*P_hat;
% Proj=PQP*PQP';
% toc

% tic
% PQPT=P_hat'*kron(P_tru',Q_tru');
% Proj=PQPT'*PQPT;
% toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% grad test
% fprintf("------------------ grad test ------------\n");
% tic
% cput=cputime;
Pro=proj_supp(UVkr);
PMA=Pro*MA;
gra=grad_opt(PMA,UVkr,U,V,MA);
% fprintf("------------------ ori method costs is %3.20f -------------\n",cputime-cput);
% toc

% tic
% cput=cputime;
% UVi=Q_hat'*inv(Sigma_hat)*P_hat'*PQkr';
% gra_opt=grad_opt_opt(R,UVi,U,V,MA);
% fprintf("------------------ ori method costs is %3.20f -------------\n",cputime-cput);
% toc

% grad_res=gra-gra_opt;
% fprintf("----------------- the residual of the grad is %3.40f----------------\n",norm(grad_res(:)));

tic
cput=cputime;
UVi=Q_hat'*inv(Sigma_hat)*P_hat'*PQkr';
PUVi=PQkr*P_hat*inv(Sigma_hat)*Q_hat;
gra_opt_opt=grad_opt_opt_opt(R,UVi,PUVi,U,V,MA);
fprintf("------------------ ori method costs is %3.20f -------------\n",cputime-cput);
toc

grad_res2=gra_opt_opt-gra;
fprintf("----------------- the residual of the grad is %3.40f----------------\n",norm(grad_res2(:)));