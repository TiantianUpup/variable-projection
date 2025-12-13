% this function calculates the value of the objective function
function eqs = fval(U,V,MA)
    % fprintf("this method calculate fval\n");
    % kr=khatrirao(U,V);
    % [mn,r]=size(kr);
    % [Q,Sigma,S] = qr(kr);
    % %[Q,Sigma,S]=svd(kr);

    % k=rank(Sigma); % the rank of the U\odot V
   
    % Ib=eye(mn-k);
    % temp=zeros(mn,mn);
    % temp(k+1:mn,k+1:mn)=Ib;
 
    % %Qtr=Q(:,k+1:mn);
    % %res=Qtr'*MA;
    P=proj_supp(U,V);
    res=P*MA;
    eqs=norm(res(:))^2;
end
