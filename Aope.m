%% this function calculates L(J)=A_1^*A_1(J)+A_2^*A_2(J)+A_3^*A_3(J)+gamma*A_4^*A_4(J)+2*lambda*A_5^*A_5(J)
function X=Aope(Aparas,paras,J,gamma,lambda)
    A=Aparas.A;
    C=Aparas.C;
    E=Aparas.E;
    G=Aparas.G;
    B2=Aparas.B2;
    B3=Aparas.B3;
    tildeU_1=Aparas.tildeU_1;
    tildeV_1=Aparas.tildeV_1;
    
    m=paras.m;
    n=paras.n;
    r=paras.r;

    X1=zeros(r,r);
    X2=zeros(m-r,r);
    X3=zeros(r,r);
    X4=zeros(n-r,r);

    J1=J(1:r,1:r);
    J2=J(r+1:m,1:r);
    J3=J(m+1:m+r,1:r);
    J4=J(m+r+1:m+n,1:r);
    
    S1=A'*(khatrirao(tildeU_1,J3)+khatrirao(J1,tildeV_1))*C;
    S2=khatrirao(tildeU_1,J4)*C;
    S3=khatrirao(J2,tildeV_1)*C;

    Ktemp_1=khatrirao(tildeU_1,J3)+khatrirao(J1,tildeV_1);
    Ktemp_2=khatrirao(tildeU_1,J4);
    Ktemp_3=khatrirao(J2,tildeV_1);
    %K=E*Ktemp_1'*G+E*Ktemp_2'*B2+E*Ktemp_3'*B3;
    K=E*(Ktemp_1'*G+Ktemp_2'*B2+Ktemp_3'*B3);

    Temp1=A*S1*C'+gamma*G*K'*E;
    Temp2=S3*C'+gamma*B3*K'*E;
    Temp3=S2*C'+gamma*B2*K'*E;

    Ir=eye(r);
    Im=eye(m-r);
    In=eye(n-r);

    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);
        X1(:,i)=kron(Ir,tildeV_1_i')*Temp1(:,i);
        X2(:,i)=kron(Im,tildeV_1_i')*Temp2(:,i);
        X3(:,i)=kron(tildeU_1_i',Ir)*Temp1(:,i);
        X4(:,i)=kron(tildeU_1_i',In)*Temp3(:,i);
    end    

    X=[X1+2*lambda*J1;X2+2*lambda*J2;X3+2*lambda*J3;X4+2*lambda*J4];
end    