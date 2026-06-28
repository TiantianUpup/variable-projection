%% this function calculates L(J)=A_1^*A_1(J)+A_2^*A_2(J)+A_3^*A_3(J)+gamma*A_4^*A_4(J)+2*mu*A_5^*A_5(J)
function X=Aope_simple(Aparas,paras,J,gamma)
    M=Aparas.M;
    MAC=Aparas.MAC;
    W2W2T=Aparas.W2W2T;
    
    tildeU_1=Aparas.tildeU_1;
    tildeV_1=Aparas.tildeV_1;

    m=paras.m;
    n=paras.n;
    r=paras.r;
   
    J1=J(1:r,1:r);
    J2=J(r+1:m,1:r);
    J3=J(m+1:m+r,1:r);
    J4=J(m+r+1:m+n,1:r);

    CCT=M*MAC;

    X1=zeros(r,r);
    X2=zeros(m-r,r);
    X3=zeros(r,r);
    X4=zeros(n-r,r);

    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);

        CCi=CCT(:,i);
        
        % CCi_diag=diag(CCi);
        % Temp_i=W2W2T*reshape(J3*CCi_diag*tildeU_1'+tildeV_1*CCi_diag*J1',[r^2,1]);
        
        Temp_i=W2W2T*reshape(J3.*CCi'*tildeU_1'+tildeV_1.*CCi'*J1',[r^2,1]);
        
        Temp_mat=reshape(Temp_i,[r,r]);
       
        X1(:,i)=Temp_mat'*tildeV_1_i;  
        X2(:,i)=J2*(CCi.*(tildeV_1'*tildeV_1_i));
        X3(:,i)=Temp_mat*tildeU_1_i;
        X4(:,i)=J4*(CCi.*(tildeU_1'*tildeU_1_i)); 
    end    
    
    X=[X1;X2;X3;X4];

    if gamma~=0
        R1=zeros(r,r);
        R2=zeros(m-r,r);
        R3=zeros(r,r);
        R4=zeros(n-r,r);

        B1=Aparas.B1;
        B3=Aparas.B3;
        B2=Aparas.B2;
        E=Aparas.E;
        W2=Aparas.W2;

        G=W2*B1;
        Ktemp_1=khatrirao(tildeU_1,J3)+khatrirao(J1,tildeV_1);
        Ktemp_2=khatrirao(tildeU_1,J4);
        Ktemp_3=khatrirao(J2,tildeV_1);
        %K=E*Ktemp_1'*G+E*Ktemp_2'*B2+E*Ktemp_3'*B3;
        K=E*(Ktemp_1'*G+Ktemp_2'*B2+Ktemp_3'*B3);
      
        GKTE=G*K'*E;
        B3KTE=B3*K'*E;
        B2KTE=B2*K'*E;

        for i=1:r
            tildeU_1_i=tildeU_1(:,i);
            tildeV_1_i=tildeV_1(:,i);

            Temp1_mat=reshape(GKTE(:,i),[r,r]);
            Temp2_mat=reshape(B3KTE(:,i),[r,m-r]);

            R1(:,i)=Temp1_mat'*tildeV_1_i;
            R2(:,i)=Temp2_mat'*tildeV_1_i;
            R3(:,i)=Temp1_mat*tildeU_1_i;
            R4(:,i)=reshape(B2KTE(:,i),[n-r,r])*tildeU_1_i;
        end    

        R=[R1;R2;R3;R4];
        X=X+gamma*R;
    end
end    