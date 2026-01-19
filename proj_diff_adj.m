function X=proj_diff_adj(U,V,H)
    [m,r]=size(U);
    [n,r]=size(V);

    Im=eye(m);
    In=eye(n);

    X1=zeros(m,r);
    X2=zeros(n,r);

    UVi=sym_inv(khatrirao(U,V));
    Z=proj_supp(U,V)*(H+H')*UVi';

    for i=1:r
        Zi=Z(:,i);
        Ui=U(:,i);
        Vi=V(:,i);
        X1(:,i)=kron(Im,Vi')*Zi;
        X2(:,i)=kron(Ui',In)*Zi;
    end  

    X=[X1;X2];
end    
