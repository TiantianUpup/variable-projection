addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")

% disp(rand(3,3));


% A=[2,2;
%    2,2];
% B=[3,3;
%    3,3;
%    3,3];
% C=[A;B];
% disp(C);   

% A=[1,2;
%    2,1;
%    1,3];
% B=[4,1,6;
%    2,1,1];
% C=A*B;
% disp(C);
U=[1,1;
   1,2;
   1,3];
V=[2,3;
   2,2;
   2,1];

W=[3,1;
   2,3;
   4,4];

MA=khatrirao(U,V)*W';
fprintf("the size of MA is\n"); 
disp(size(MA));
A=khatrirao(U,V);   
     
[Q,Sigma,S]=qr(A);
S=S';
r=rank(A);
[m,n]=size(A);

Tau=zeros(n,m);
A11=Sigma(1:r,1:r);
Tau(1:r,1:r)=inv(A11);
Ai=S'*Tau*Q';

disp(S);
% disp(A);
% disp(A*Ai*A);

res=A-A*Ai*A;
disp(norm(res(:)));



% AAi=A*Ai;
% disp(AAi');

% disp(AAi);

% res=AAi-AAi';
% fprintf("res is %d\n",norm(res(:)));

% I=eye(m,m);
% P=A*Ai;
% Pc=I-P;
% Pcf=Pc*MA;
% fprintf("the size of PCf is\n"); 
% disp(size(Pcf));

% fprintf("value is %3.4f\n",norm(Pcf(:)));

%Pci=proj_supp(U,V);
%disp(Pci);


%H=Pc+Pci;

%disp(H);
% Pcif=Pci*MA;
% disp(size(Pcif));
% disp(Pcif);
% fprintf("value is %3.4f\n",norm(Pcf(:)));
% res=Pc-Pci;
% fprintf("residual is %3.12f\n",norm(res(:)));

% pres=Pcf-Pcif;
% fprintf("residual is %3.12f\n",norm(pres(:)));

