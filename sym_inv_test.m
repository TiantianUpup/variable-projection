m=10;
n=8;
r=3;

U=rand(m,r);
V=rand(n,r);
addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")
UVkr=khatrirao(U,V);
UVkri=sym_inv(UVkr);
I=eye(m*n);
Pc=I-UVkr*UVkri;
Pc_2=proj_supp(U,V);
res=Pc-Pc_2;
fprintf("residual is %3.20f\n",norm(res(:)));

