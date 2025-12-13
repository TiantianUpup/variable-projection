% initial point
x_1=12;
x_2=15;
x=[x_1;x_2];

g_temp_1=x_1^2+x_2^2-1;
g_temp_2=x_1+x_2-2;
gval=g_temp_1^2+g_temp_2^2;

% solve the normal equation J'Jd=-J'f by the cg method
% A=J'Jd, B=-J'f
A=[4*x_1^2+1,4*x_1*x_2+1;
   4*x_1*x_2+1,4*x_2^2+1];
b=[2*x_1^3+2*x_1*x_2^2-x_1+x_2-2;
   2*x_1^2*x_2+2*x_2^3+x_1-x_2-2]; 
grad=b;     
b=-b;

gn_tol=1e-4;
gn_itmax=150;
iter=0;

paras.itmax=50;
paras.tol=1e-6;

while (norm(grad)>gn_tol && iter<gn_itmax)
   iter=iter+1;
   d=gn_cg(A,b,paras);
   
   inn=grad'*d;

   fprintf("inn is %3.4f\n",inn);

   % new iteration point
   x=x+d;

   % update the gval
   x_1=x(1);
   x_2=x(2);
   g_temp_1=x_1^2+x_2^2-1;
   g_temp_2=x_1+x_2-2;
   gval=g_temp_1^2+g_temp_2^2;
   fprintf("iter is %d, gval_cur is %3.4f\n",iter, gval);

   % update A and b, and gradient
   A=[4*x_1^2+1,4*x_1*x_2+1;
   4*x_1*x_2+1,4*x_2^2+1];
   b=[2*x_1^3+2*x_1*x_2^2-x_1+x_2-2;
   2*x_1^2*x_2+2*x_2^3+x_1-x_2-2]; 
   grad=b;     
   b=-b; 

   fprintf("the norm of the gradient is %3.4f\n",norm(grad));
end

disp(x);
