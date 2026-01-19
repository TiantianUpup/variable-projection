function eqs = fval_tilde(P1,Q1,W1,MA)
    %P=proj_supp(U,V);
    % A=khatrirao(U,V);
    % P=proj_supp(A);
    % res=P*MA;
    % eqs=0.5*norm(res(:))^2;
    %P=[p_1,...,p_r] vec(A*P) 

    [m,r]=size(P1);
    [n,r]=size(Q1);

    T=zeros(m*n,r);
    for i=1:r
        wi=W1(:,i);
        % mat(wi)=reshape(w_i,[r,r]);
        % Temp_i=Q1*mat(wi)*P1';
        % T_i=vec(Temp_i)=reshape(Temp_i,[m*n,1])
        T(:,i)=reshape(Q1*reshape(wi,[r,r])*P1',[m*n,1]);

    end
    MAT=MA'*T;    
    eqs=-0.5*sum(MAT(:).^2);
end