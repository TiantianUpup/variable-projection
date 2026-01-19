%%% this function calculates the projection supplement:P_{A}^{\prep}
function P = proj_supp(A)
    % kr=khatrirao(U,V);
    % [mn,r]=size(kr);
    
    [mn,r]=size(A);
    [Q,Sigma,S] = qr(A);
    %[Q,Sigma,S]=svd(kr);

    k=rank(Sigma); % the rank of the U\odot V

    Q1=Q(1:mn,1:k);
    I=eye(mn);
    P=I-Q1*Q1';
end