%% this function calculates L(J)=A_1^*A_1(J)+A_2^*A_2(J)+A_3^*A_3(J)+gamma*A_4^*A_4(J)+2*mu*A_5^*A_5(J)
function X=Aope_opt(Aparas,paras,J,gamma,mu,lambda)
    A=Aparas.A;
    C=Aparas.C;
    E=Aparas.E;
    G=Aparas.G;
    B2=Aparas.B2;
    B3=Aparas.B3;
    % if isfield(Aparas, {'AS','P','Q'})
    %     AS=Aparas.AS;
    %     P=Aparas.P;
    %     Q=Aparas.Q;
    % end

    tildeU_1=Aparas.tildeU_1;
    tildeV_1=Aparas.tildeV_1;
 
    m=paras.m;
    n=paras.n;
    r=paras.r;

    J1=J(1:r,1:r);
    J2=J(r+1:m,1:r);
    J3=J(m+1:m+r,1:r);
    J4=J(m+r+1:m+n,1:r);

    % These codes can be optimized, khatrirao(tildeU_1,J4)=T1, khatrirao(J2,tildeV_1)=T2 
    U1_J3_kr=khatrirao(tildeU_1,J3);
    J1_V1_kr=khatrirao(J1,tildeV_1);
    U1_J4_kr=khatrirao(tildeU_1,J4);  % Ktemp_2
    J2_V1_kr=khatrirao(J2,tildeV_1);  % Ktemp_3

    S1=A'*(U1_J3_kr+J1_V1_kr)*C;
    S2=U1_J4_kr*C;
    S3=J2_V1_kr*C;

    Ktemp_1=U1_J3_kr+J1_V1_kr;
   
    % tic
    % K_opt=E*Ktemp_1'*G+E*U1_J4_kr'*B2+E*J2_V1_kr'*B3;
    % toc

    K=E*(Ktemp_1'*G+U1_J4_kr'*B2+J2_V1_kr'*B3);
    
    Temp1=A*S1*C'+gamma*G*K'*E;
    Temp2=S3*C'+gamma*B3*K'*E;
    Temp3=S2*C'+gamma*B2*K'*E;

    X1=zeros(r,r);
    X2=zeros(m-r,r);
    X3=zeros(r,r);
    X4=zeros(n-r,r);
   
    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);

        Temp1_mat=reshape(Temp1(:,i),[r,r]);
        Temp2_mat=reshape(Temp2(:,i),[r,m-r]);
        
        X1(:,i)=Temp1_mat'*tildeV_1_i;  
        X2(:,i)=Temp2_mat'*tildeV_1_i;
        X3(:,i)=reshape(Temp1(:,i),[r,r])*tildeU_1_i;
        X4(:,i)=reshape(Temp3(:,i),[n-r,r])*tildeU_1_i; 
    end    
    
    X=[X1+mu*J1;X2+mu*J2;X3+mu*J3;X4+mu*J4];
 

    % X=[X1+2*mu*J1;X2+2*mu*J2;X3+2*mu*J3;X4+2*mu*J4];

    % X=[X1+mu*J1;X2+mu*J2;X3+mu*J3;X4+mu*J4];

    % X_res=X_opt-X;
    % fprintf("X residual is %3.20f\n",norm(X_res(:)));
    % if lambda

    %     S11=S(1:m,1:m);
    %     S12=S(1:m,m+1:m+n);
    %     S21=S(m+1:m+n,1:m);
    %     S22=S(m+1:m+n,m+1:m+n);

    %     Temp=[P'*S11*P*J12+P'*S12*Q*J34;
    %           Q'*S21*P*J12+Q'*S22*Q*J34];

    %     X=X+Temp;      
    % end    

    %X=X_opt;
end    