clc
clear all
close all

addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")

% m=6;
% n=5;
% r=3;

m=160;
n=200;
r=68;

% m=10;
% n=12;
% r=5;


H=permutation(m,n,r);
%disp(H);
H1 = H(1:m*n,1:r^2);
H21 = H(1:m*n,(r^2+1):n*r);
H22 = H(1:m*n,n*r+1:(n+m)*r-r^2);
H23 = H(1:m*n,(n+m)*r-r^2+1:m*n);

%%%% (m,n,r)=(10 12 5)
% U=[1,2,1,1,5;
%    2,5,2,3,5;
%    3,4,5,6,5;
%    4,7,8,9,5;
%    5,8,7,3,5;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0];

% V=[4,3,2,1,5;
%    3,2,2,4,5;
%    2,4,6,8,5;
%    1,0,2,4,5;
%    7,9,8,8,5; 
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0;
%    0,0,0,0,0];   

% %%% (m,n,r)=(6 5 3)
% U=[1,2,1;
%    2,5,2;
%    3,4,5;
%    0,0,0;
%    0,0,0;
%    0,0,0];

% V=[4,3,2;
%    3,2,2;
%    2,4,6;
%    0,0,0;
%    0,0,0];   
   
% Utilde=U(1:r,1:r);
% Vtilde=V(1:r,1:r);


% %%%% (m,n,r)=(10 8 4)
% X=[1,4,2,3;
%    2,8,5,4;
%    3,2,11,13;
%    4,16,3,7;
%    5,6,7,8;
%    6,9,7,5;
%    7,9,11,12;
%    8,9,10,11;
%    9,4,6,7;
%    10,9,8,7];
  
% Y=[1,3,5,7;
%    2,7,4,1;
%    3,1,2,2;
%    4,3,2,6;
%    5,4,3,2;
%    6,4,2,8;
%    7,3,4,5;
%    8,2,2,2];  



%%% (m,n,r)=(6 5 3)
U=zeros(m,r);
Utilde=rand(r,r);
U(1:r,1:r)=Utilde;
V=zeros(n,r);
Vtilde=rand(r,r);
V(1:r,1:r)=Vtilde;

X=rand(m,r);
Y=rand(n,r);  

%%%% (m,n,r)=(10 12 5)
% X=[1,4,2,3,4;
%    2,8,5,4,6;
%    3,2,11,7,8;
%    4,16,3,9,7;
%    5,6,7,7,8;
%    6,9,7,9,9;
%    7,5,7,3,1;
%    8,5,5,3,14;
%    9,4,2,2,5;
%    10,6,7,8,9];
  
% Y=[1,3,5,7,9;
%    2,7,4,6,7;
%    3,1,2,7,8;
%    4,3,2,1,4;
%    5,4,3,2,1;
%    6,5,4,3,2;
%    7,2,1,6,9;
%    8,9,3,5,7;
%    9,4,3,1,2;
%    10,9,8,7,6;
%    11,3,2,7,8;
%    12,8,7,6,5]; 


X1=X(1:r,1:r);
X2=X(r+1:m,1:r);    
Y1=Y(1:r,1:r);
Y2=Y(r+1:n,1:r);

% fprintf("============= H1 left ==============\n");
% disp(H1'*khatrirao(U,Y));
% fprintf("============= H1 right ==============\n");
% disp(khatrirao(Utilde,Y1));
R1=H1'*khatrirao(U,Y)-khatrirao(Utilde,Y1);
fprintf("residual is %3.8f\n",norm(R1(:)));

% fprintf("============= H21 left ==============\n");
% disp(H21'*khatrirao(U,Y));
% fprintf("============= H21 right ==============\n");
% disp(khatrirao(Utilde,Y2));
R2=H21'*khatrirao(U,Y)-khatrirao(Utilde,Y2);
fprintf("residual is %3.8f\n",norm(R2(:)));

% fprintf("============= H22 left ==============\n");
% disp(H22'*khatrirao(U,Y));
R3=H22'*khatrirao(U,Y);
fprintf("residual is %3.8f\n",norm(R3(:)));

% fprintf("============= H23 left ==============\n");
% disp(H23'*khatrirao(U,Y));
R4=H23'*khatrirao(U,Y);
fprintf("residual is %3.8f\n",norm(R4(:)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf("============= H1 left ==============\n");
% disp(H1'*khatrirao(X,V));
% fprintf("============= H1 right ==============\n");
% disp(khatrirao(X1,Vtilde));
R1=H1'*khatrirao(X,V)-khatrirao(X1,Vtilde);
fprintf("residual is %3.8f\n",norm(R1(:)));

% fprintf("============= H21 left ==============\n");
% disp(H21'*khatrirao(X,V));
R2=H21'*khatrirao(X,V);
fprintf("residual is %3.8f\n",norm(R2(:)));


% fprintf("============= H22 left ==============\n");
% disp(H22'*khatrirao(X,V));
% fprintf("============= H22 right ==============\n");
% disp(khatrirao(X2,Vtilde));
R3=H22'*khatrirao(X,V)-khatrirao(X2,Vtilde);
fprintf("residual is %3.8f\n",norm(R3(:)));

% fprintf("============= H23 left ==============\n");
% disp(H23'*khatrirao(X,V));
R4=H23'*khatrirao(X,V);
fprintf("residual is %3.8f\n",norm(R4(:)));
