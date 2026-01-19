function eqs = fval_tilde_opt(WTW,P1,Q1,U,V,W1)
    fprintf("-------------- this is fval_tilde_opt method -------------------\n");
    [k,~]=size(WTW);
    [~,r]=size(P1);
    W=cell(r);

    for i=1:r
        % wi is the i-th column of W
        wi=W1(:,i);
        w_mat=reshape(wi,[r,r]);
        W{i}=w_mat;
    end

    N=zeros(r,k);

    % fprintf("------------------ method 1 --------------------\n");
    % tic
    % for i=1:r
    %     wi=W1(:,i);
    %     W{i}=Q1*reshape(wi,[r,r])*P1';
    %     for j=1:k
    %         vj=V(:,j);
    %         uj=U(:,j);
    %         N(i,j)=vj'*W{i}*uj;
    %     end    
    % end
    % toc
    % fprintf("------------------ method 1 end --------------------\n");

    % N2=zeros(r,k);
    % fprintf("------------------ method 2 --------------------\n");
    % tic
    for i=1:r
        wi=W1(:,i);
        W{i}=reshape(wi,[r,r]);
        for j=1:k
            vj=V(:,j);
            uj=U(:,j);
            N(i,j)=(vj'*Q1)*W{i}*(P1'*uj);
        end    
    end
    % toc
    % fprintf("------------------ method 2 end --------------------\n");

    %N2TN2=N2'*N2; 
    NTN=N'*N;
    
    %N_res=NTN-N2TN2;

    %fprintf("----------------------- the residual of N is %3.20f -------------\n",norm(N_res));
    eqs=-0.5*sum(WTW(:).*NTN(:));
end