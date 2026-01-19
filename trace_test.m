% A=[1,2;
%    3,3;
%    11,201;
%    2,34];
% B=[3,4;
%    25,56;
%    98,44;
%    45,6];
   
% disp(trace(A'*B));
% disp(trace(A*B'));   

% A=[1,2;
%    3,4];
% disp(A);
% A=A';
% disp(A);   


% A=rand(4,3);
% [P1, Sigma1, S1] = qr(U);
% fprintf("first\n");
% disp(P1);
% disp(Sigma1);
% fprintf("second\n");
% [P2, Sigma2, S2] = qr(U);
% disp(P2);
% disp(Sigma2);

% res1=P1-P2;
% res2=Sigma1-Sigma2;
% fprintf("res is %3.20f\n",norm(res1(:)));
% fprintf("res is %3.20f\n",norm(res2(:)));

% for i=0:1
%    A=[i,i;
%       i,i];

%    B=[i+1,i+1;
%       i+1,i+1];
%    C=[A;B];   
%    fprintf("the norm of C is %3.8f\n",norm(C(:))^2);      
% end   

A=[1,4;
   2,5;
   3,6];
disp(A(:));  

a=A(:);
A=reshape(a,3,2);
disp(A);