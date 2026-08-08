%% this function calculates A^*A(C)
function X = AadjointA (Aparas,paras,J)                         
    W1=Aparas.W1;
    CCT=Aparas.CCT;

    tildeU_1=Aparas.tildeU_1;
    tildeV_1=Aparas.tildeV_1;

    r=paras.r;
   
    J1=J(1:r,1:r);
    J3=J(r+1:2*r,1:r);

    X1=zeros(r,r);
    X3=zeros(r,r);
   
    for i=1:r
        CCi=CCT(:,i);
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);
        a = reshape(J3.*CCi'*tildeU_1'+tildeV_1.*CCi'*J1',[r*r,1]);
        a_mat=reshape(a-W1*(W1'*a),[r,r]);
        
        X1(:,i)=a_mat'*tildeV_1_i;  
        X3(:,i)=a_mat*tildeU_1_i;
    end    

    %X_opt=[X1;X2;X3;X4];  
    X=[X1;X3];   
end     