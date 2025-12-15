function X=Lope(Lparas,paras,C)
    U=Lparas.U;
    V=Lparas.V;
    MA=Lparas.MA;
    
    m=paras.m;
    n=paras.n;
    r=paras.r;

    C1=C(1:m,1:r);
    C2=C(m+1:m+n,1:r);

    X1=zeros(m,r);
    X2=zeros(n,r);

    kri=sym_inv(khatrirao(U,V));
    P=proj_supp(U,V);
    Htemp=P*(khatrirao(U,C2)+khatrirao(C1,V))*kri;
    H=-(Htemp+Htemp')*MA;

    Temp=-H*MA';
    Z=P*(Temp+Temp')*kri';
    
    Im=eye(m);
    In=eye(n);

    for i=1:r
        Ui=U(:,i);
        Vi=V(:,i);
        Zi=Z(:,i);
        X1(:,i)=kron(Im,Vi')*Zi;
        X2(:,i)=kron(Ui',In)*Zi;
    end

    X=[X1;X2];
end    