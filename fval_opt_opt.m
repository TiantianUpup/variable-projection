function eqs = fval_opt_opt(PMA)
    fprintf("------------ this is fval_opt method --------------\n");
    %P=proj_supp(U,V);
    % A=khatrirao(U,V);
    % P=proj_supp(A);
    % res=P*MA;
    % eqs=0.5*norm(res(:))^2;

    eqs=0.5*norm(PMA(:))^2;
end