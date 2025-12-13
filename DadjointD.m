%% this function calculates A^*A(C)
function X = DadjointD (tildeU_1,tildeV_1,MMA,WL,DF,H,paras, C)                         
    r=paras.r;
    m=paras.m;
    n=paras.n;
    X=zeros(m+n,r);
    X1=zeros(r,r);
    X2=zeros(m-r,r);
    Y1=zeros(r,r);
    Y2=zeros(n-r,r);

    I1=eye(r,r);
    I2=eye(m-r,m-r);
    I3=eye(n-r,n-r);
 
    C1=C(1:r,1:r);
    C2=C(r+1:m,1:r);
    C3=C(m+1:m+r,1:r);
    C4=C(m+r+1:m+n,1:r);

    Temp=B'*(khatrirao(tildeU_1,C3)+khatrirao(C1,tildeV_1))+C'*khatrirao(tildeU_1,Y1)+D'*khatrirao(C2,tildeV_1);
    Temp1=WL'*WL*(khatrirao(tildeU_1,C3)+khatrirao(C1,tildeV_1))*MMA*MMA'+B*Temp*A'*A;
    Temp2=khatrirao(C2,tildeV_1)*MMA*MMA'+D*Temp*A'*A;
    Temp3=khatrirao(tildeU_1,C4)*MMA*MMA'+C*Temp*A'*A;

    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);
        X1(:,i)=kron(I1,tildeV_1_i')*Temp1(:,i);
        X2(:,i)=kron(I2,tildeV_1_i')*Temp2(:,i);
        Y1(:,i)=kron(tildeU_1_i',I1)*Temp1(:,i);
        Y2(:,i)=kron(tildeU_1_i',I3)*Temp3(:,i);
    end 

    X(1:r,1:r)=X1;
    X(r+1:m,1:r)=X2;
    X(m+1:m+r,1:r)=Y1;
    X(m+r+1:m+n,1:r)=Y2;
end     