% this function calculates the value of the objective function
function eqs = fval(U,V,MA)
    %P=proj_supp(U,V);
    A=khatrirao(U,V);
    P=proj_supp(A);
    res=P*MA;
    
    eqs=0.5*norm(res(:))^2;
end
