%% this function calculates A^*A(C)
function X = AadjointA (tildeU_1,tildeV_1,MMA,WL,paras, C)                         
    r=paras.r;
    X = zeros(2*r,r);
    X1=zeros(r,r);
    X2=zeros(r,r);
    I = eye(r);

    % tildeU_1 = operator_paras.tildeU_1;
    % tildeV_1 = operator_paras.tildeV_1;
    % MMA = operator_paras.MMA;
    % WL = operator_paras.WL; % WL: W*L

    C1=C(1:r,1:r);
    C2=C(r+1:2*r,1:r);

    temp=WL*WL'*(khatrirao(tildeU_1,C2)+khatrirao(C1,tildeV_1))*(MMA*MMA');
    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);
        X1(:,i)=kron(I,tildeV_1_i')*temp(:,i);
        X2(:,i)=kron(tildeU_1_i',I)*temp(:,i);
    end 

    X(1:r,1:r)=X1;
    X(r+1:2*r,1:r)=X2;
end     