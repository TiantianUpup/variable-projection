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
m=30;
n=40;
r=20;
% U=[1,2;
%    3,4;
%    5,6];
% Y=[2,3;
%    7,9];
% Z=[2,4;
%    6,8;
%    10,12;
%    14,16;
%    18,20;
%    1,5];

% A=[1,2,3,4,5,6;
%    2,4,6,8,10,12;
%    3,1,2,5,4,0;
%    4,0,2,7,2,4;
%    5,7,9,8,2,3;
%    6,5,9,8,7,6];   
% B=[6,9;
%    5,3];
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
% L=A*khatrirao(U,Y)*B;
% val_1=trace(L'*Z);
% fprintf("val_1 is %d\n",val_1);
% tildeZ=A'*Z*B';
% Temp=zeros(n,r);
% I=eye(n);
% for i=1:r
%     Ui=U(:,i);
%     tildeZi=tildeZ(:,i);
%     Temp(:,i)=kron(Ui',I)*tildeZi;
% end    
% val_2=trace(Y'*Temp);   
% fprintf("val_2 is %d\n",val_2);
% U=rand(m,r);
% V=rand(n,r);
% R1=khatrirao(U,V);
% R2=kr(U,V);


% a=[1;2;3;4;5;6;7;8;9];
% a_mat=reshape(a,[3,3]);
% disp(a_mat);

% a = 0.1 + 0.2 - 0.3;
% disp(a);


% r = 3;
% m = 5;
% i = 1;

% % 生成随机但可重复的数据
% rng(1);
% Temp1 = randn(r^2, 10);
% tildeV_1 = randn(r, 10);
% Ir = eye(r);
% Im = eye(m);

% % 你的原始计算
% tildeV_1_i = tildeV_1(:,i);
% Temp1_mat = reshape(Temp1(:,i), [r, r]);

% X1_opt_i = Temp1_mat' * tildeV_1_i;
% X1_i = kron(Ir, tildeV_1_i') * Temp1(:,i);

% residual = X1_opt_i - X1_i;
% disp('Residual:');
% disp(residual);
% disp(['Norm of residual: ', num2str(norm(residual))]);
% disp(['Relative norm: ', num2str(norm(residual)/norm(X1_opt_i))]);


% X = ktensor({rand(4,2),rand(2,2),rand(3,2)}) %<-- Data.
% Y = ktensor({rand(4,2),rand(2,2),rand(3,2)}) %<-- More data.

% disp("Z is");
% Z=tensor(X)-tensor(Y)

m=100; n=150; r=35;
% X1=[1,2;
%     3,4];
% X2=[5,6];
% Y1=[7,8;
%     9,10];
% Y2=[11,12;
%     13,14];

% X1=rand(r,r);
% X2=rand(m-r,r);
% Y1=rand(r,r);
% Y2=rand(n-r,r);

% Z=[X1;X2;Y1;Y2];
% w=[X1(:);
%    X2(:);
%    Y1(:);
%    Y2(:)];
% z=Z(:);
% P=perm_cg(m,n,r);
% Pz=P*z;

% res=Pz-w;

% disp(Pz);
% disp("----------------------------------------------------------");
% disp(w);
%fprintf("---------------------- the residual is %3.50f -----------------------\n",norm(res));


% A=[1,2,3;
%    4,2,6;
%    7,8,4];
% disp(1./diag(A));   

% n = 50; 
% A = randn(n); 
% A = A' * A; 
% A = A + eye(n) * 0.1; 
% b=A*ones(n,1);

% cg_paras.itmax=1500;
% cg_paras.tol=1e-6;
% cg_paras.lambda=0;
% runhist=pcg(A,b,cg_paras);

% runhist_mcg=mcg(A,b,cg_paras);
% x=runhist.x;
% disp(x);

% m=2;
% n=2;
% r=2;

% P=rand(r,(m+n)*r);
% v=rand(n,1);
% result_1=kron(P,v);

% result_2=zeros(n*r,(m+n)*r);
% for i=1:(m+n)*r
%     result_2(:,i)=kron(P(:,i),v);
% end

% res=result_1-result_2;
% fprintf("the residual is %3.80f\n",norm(res(:)));

% r=3;
% W=cell(r);

% for i=1:r
%     W{i}=i*ones(r);
% end    

% fprintf("--------------------- W{1} is ------------------------\n");
% disp(W{1});
% fprintf("--------------------- W{2} is ------------------------\n");
% disp(W{2});
% fprintf("--------------------- W{3} is ------------------------\n");
% disp(W{3});

% W=[1,2;
%    3,4;
%    5,6];
% disp(W(1,:));
% disp(W(2,:));
% disp(W(3,:));  

% m=3;
% n=4;
% r=2;
% X=rand(m,r);
% Y=rand(n,r);
% v=rand(r,1);

% result_1=khatrirao(X,Y)*v;
% Temp=Y*diag(v)*X';
% result_2=Temp(:);
% disp(result_1);
% disp(result_2);

% m=3;
% r=4;
% s=2;
% t=4;
% X=rand(m,r);
% Y=rand(r,s);
% A=rand(n,t);
% B=rand(t,s);

% result_1=khatrirao(X*Y,A*B);
% result_2=kron(X,A)*khatrirao(Y,B);
% res=result_1-result_2;
% fprintf("------------------- the residual is ------------------%3.20f\n",norm(res(:)));

% x=0|5;
% disp(mod(1,5));
% disp(1|5);
% disp(2|5);
% disp(3|5);
% disp(4|5);
% disp(5|5);
% disp(6|5);
% r=5;
% v=rand(r,1);
% R=zeros(r,r);
% A1=rand(r,r);
% A2=rand(r,r);
% A3=rand(r,r);
% A4=rand(r,r);
% A5=rand(r,r);
% A=[A1;A2;A3;A4;A5];
% I=eye(r);
% result_1=kron(I,v)*A;
% % for i=1:r
% %     R(i,1:end)=;
% % end
% Arow=[A1,A2,A3,A4,A5];
% R_temp=v'*Arow;
% result_2=reshape(R_temp,[r,r]);
% res=result_1-result_2;
% fprintf("the residual is %3.4f\n",res);


% A=[1,2;
%    3,4];
% a=diag(A);   
% disp(a);

% disp(sum(A,1));

% m=4;
% n=3;
% P=rand(m,m);
% Q=rand(n,n);
% u=rand(m,1);
% v=rand(n,1);
% alpha=rand(m*n,1);
% uTP=u'*P;
% vTQ=v'*Q;
% res_1=kron(uTP,vTQ)*alpha;
% res_2=vTQ*reshape(alpha,[n,m])*uTP';
% res=res_1-res_2;
% fprintf("----------- the residual of result is %3.20f ----------------\n",res);

% m=30;
% n=40;
% r=2;
% X=rand(n,r);
% Y=rand(m,r);
% M=rand(n,m);

% res_1=zeros(r,1);
% tic
% for i=1:r
%    Xi=X(:,i);
%    Yi=Y(:,i);
%    res_1(i)=Xi'*M*Yi;
% end   
% toc

% tic
% res_2= sum((X' * M) .* Y', 2);
% toc

% res=res_1-res_2;

% fprintf("------------------- the residual is %3.20f ---------------\n",norm(res));


% A=[
%    1,2,3;
%    -0.1,-0.2,1;
%    0,1,1;
% ];
% count = sum(all(A >= 0, 2));
% fprintf("-------- count is %3.2f ----------\n", count);
% fprintf("the max nonnegative entry of A=%3.4f\n",max(A(A >= 0)));

% a=[1,2,3,4];
% b=[4,3,2,1];
% disp(khatrirao(a,b));

% a=[1,2,3,4,5,6,7,8,9];
% A=reshape(a,[3,3]);
% disp(A');

a=3.1234567;
b=round(a,2);
disp(b);