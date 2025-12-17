% this function calculates the value of the objective function
function eqs = fval_2(U,V,P,Q,F,Di,HWJ,MA)
    [m,r]=size(U);
    [n,r]=size(V);
    I=eye(m*n);
   
    % UVkr=khatrirao(U,V);
    % Temp=kron(P,Q)*HWJ;
    % UVkri=F'*Di*Temp';

    % Pc=I-UVkr*UVkri;
    % res=Pc*MA;
    
    
    PQkr=kron(P,Q);
    Pc=PQkr*(I-HWJ*HWJ')*PQkr';

    fprintf("the norm of the Pc is %3.8f\n",norm(Pc(:)));
    Pc_2=proj_supp(U,V);
    fprintf("the norm of the Pc_2 is %3.8f\n",norm(Pc_2(:)));
    res=Pc*MA;
    eqs=norm(res(:))^2;
end
