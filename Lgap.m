function eqs = Lgap(U,V,A,B,MA)
    Pc=proj_supp(U,V);

    % dP_(U\odot V)(H_1,H_2)
    UVkr=khatrirao(U,V);
    dPtemp=Pc*(khatrirao(U,B)+khatrirao(A,V))*sym_inv(UVkr);
    dP=dPtemp+dPtemp';
    dPMA=dP*MA;
    eqs_1=norm(dPMA(:))^2;
    PMA=Pc*MA;
    eqs_2=trace(PMA'*dPMA);
    eqs=-0.5*eqs_1+eqs_2;
end    