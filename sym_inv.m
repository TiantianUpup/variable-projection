%%%% this function calculates the symmetric generalized inverse of a matrix
function I = sym_inv(X)
    [m,n]=size(X);
    [Q,Sigma,S]=qr(X); % X=Q*Sigma*S'
    r=rank(Sigma);
    A11=Sigma(1:r,1:r);
    Z=zeros(n,m);
    Ai=inv(A11);
    Z(1:r,1:r)=Ai;
    
    I=S*Z*Q';

    % fprintf("sys_inc cost is %3.20f\n",cputime-cput);
    
    % cput_opt=cputime;
    % tic
    % Qt=Q';
    % Q1=Qt(1:r,1:r);
    % Q2=Qt(1:r,r+1:m);
    % Sub=S(1:n,1:r)*Ai;

    % I_opt=[Sub*Q1,Sub*Q2];
    % toc
    % fprintf("sys_inc cost is %3.20f\n",cputime-cput_opt); 

    % I_res=I-I_opt;
    % fprintf("-------------------------- the residual of I is %3.20f\n",norm(I_res(:)));

end    