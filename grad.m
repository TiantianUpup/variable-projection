%%% this function calculates the gradient of f(U,V)
function G = grad(U,V,MA)
    [m,r]=size(U);
    [n,r]=size(V);

    Im=eye(m);
    In=eye(n);

    % gradient
    G1=zeros(m,r);
    G2=zeros(n,r);

    P=proj_supp(U,V);
    PMA=P*MA;
    UVi=sym_inv(khatrirao(U,V));

    % calculate Z
    Z=-(PMA*MA'+PMA*PMA')*UVi';

    for i=1:r
        Zi=Z(:,i);
        Ui=U(:,i);
        Vi=V(:,i);
        G1(:,i)=kron(Im,Vi')*Zi;
        G2(:,i)=kron(Ui',In)*Zi;
    end    

    G=[G1;G2];
end    