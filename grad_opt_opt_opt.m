%%% this function calculates the gradient of f(U,V)
%function G = grad(U,V,MA)
function G=grad_opt_opt_opt(PMA,UVi,PUVi,U,V,MA)
    fprintf("------------------------------ this is grad_opt_opt_opt method ----------------------\n");
    [m,r]=size(U);
    [n,r]=size(V);

    % Im=eye(m);
    % In=eye(n);

    % gradient
    G1=zeros(m,r);
    G2=zeros(n,r);

    % UVkr=khatrirao(U,V);
    % P=proj_supp(UVkr);

    % PMA=P*MA;
    % UVi=sym_inv(UVkr);

    % P=proj_supp(UVkr);
    % PMA=P*MA;
    % UVi=sym_inv(UVkr);

    % calculate Z
    RMA=PMA*MA';
    Z=-(RMA*UVi'+RMA*PUVi);

    % for i=1:r
    %     Zi=Z(:,i);
    %     Ui=U(:,i);
    %     Vi=V(:,i);
    %     G1(:,i)=kron(Im,Vi')*Zi;
    %     G2(:,i)=kron(Ui',In)*Zi;
    % end    

    % G=[G1;G2];

    for i=1:r
        Zi=Z(:,i);
        Ui=U(:,i);
        Vi=V(:,i);

        Z_mat=reshape(Zi,[n,m]);

        G1(:,i)=Z_mat'*Vi;
        G2(:,i)=Z_mat*Ui;
    end    

    G=[G1;G2];
end    