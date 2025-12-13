% m=100;
% n=10;
% p=10;

% r=2;
% % truth factor matrices
% U = rand(m,r);

% %disp(rank(U));
% V = rand(n,r);
% %disp(rank(V));
% W = rand(p,r);

U=[1,1;
   1,1;
   1,1];
V=[2,2;
   2,2;
   2,2];
W=[3,3;
   3,3;
   3,3];      

MA=khatrirao(U,V)*W';

fprintf("fval is %3.8f\n",fval(U,V,MA));