%%%% this function calculates the symmetric generalized inverse of a matrix
function I = sym_inv(X)
    [m,n]=size(X);
    r=rank(X);
    [Q,Sigma,S]=qr(X); % X=Q*Sigma*S'
    A11=Sigma(1:r,1:r);
    Z=zeros(n,m);
    Z(1:r,1:r)=inv(A11);
    I=S*Z*Q';
end    