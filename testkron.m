% v1 = [1; 2; 3];    % 列向量
% v2 = [4; 5];       % 列向量
% v3 = [6; 7; 8];    % 列向量

% % 连续进行 Kronecker 积
% result = kron(kron(v1, v2), v3);
% disp(result);

% A=[1,2;3,4];
% A=-A;
% disp(A);
addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
m=3;
n=2;
r=2;
U=[1,2;
   3,4;
   5,6];
Y=[2,3;
   7,9];
Z=[2,4;
   6,8;
   10,12;
   14,16;
   18,20;
   1,5];

A=[1,2,3,4,5,6;
   2,4,6,8,10,12;
   3,1,2,5,4,0;
   4,0,2,7,2,4;
   5,7,9,8,2,3;
   6,5,9,8,7,6];   
B=[6,9;
   5,3];
% val_1=trace(Z'*khatrirao(U,Y));
% fprintf("val_1 is %d\n",val_1);
% Temp=zeros(n,r);
% I=eye(n);
% for i=1:r
%     Ui=U(:,i);
%     Zi=Z(:,i);
%     Temp(:,i)=kron(Ui',I)*Zi;
% end    
% val_2=trace(Y'*Temp);   
% fprintf("val_2 is %d\n",val_2);
L=A*khatrirao(U,Y)*B;
val_1=trace(L'*Z);
fprintf("val_1 is %d\n",val_1);
tildeZ=A'*Z*B';
Temp=zeros(n,r);
I=eye(n);
for i=1:r
    Ui=U(:,i);
    tildeZi=tildeZ(:,i);
    Temp(:,i)=kron(Ui',I)*tildeZi;
end    
val_2=trace(Y'*Temp);   
fprintf("val_2 is %d\n",val_2);