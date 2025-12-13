%%% this function calculates the projection supplement:P_{U\odot V}^{\prep}
function P = proj_supp(U,V)
    kr=khatrirao(U,V);
    [mn,r]=size(kr);
    [Q,Sigma,S] = qr(kr);
    %[Q,Sigma,S]=svd(kr);

    k=rank(Sigma); % the rank of the U\odot V

    Ib=eye(mn-k);
    temp=zeros(mn,mn);
    temp(k+1:mn,k+1:mn)=Ib;

    %Qtr=Q(:,k+1:mn);
    %res=Qtr'*MA;
    P=Q*temp*Q';
end