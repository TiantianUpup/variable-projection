function X=proj_diff(U,V,A,B)
    UVkr=khatrirao(U,V);
    
    Temp=proj_supp(UVkr)*(khatrirao(U,B)+khatrirao(A,V))*sym_inv(UVkr);
    X=Temp+Temp';
end    