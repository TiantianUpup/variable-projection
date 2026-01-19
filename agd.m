%%%% This is the accelerated gradient descent method
function runhist = agd(U0,V0,MA,alpha,rho,epsilon)
    fprintf("this is the agd method\n");
    timect = cputime;

    % % initial point
    % U=U0;
    % V=V0;
    % Uy=U0;
    % Vy=V0;

    % gra = grad(Uy,Vy,MA);
    % [m,r]=size(U);
    % [n,r]=size(V);
    % grad_U=gra(1:m,1:r);
    % grad_V=gra(m+1:m+n,1:r);
  
    % fun_val=fval(Uy,Vy,MA); 
    % fprintf("initial objective function value is %3.4f\n",fun_val);
    % iter=0;

    % lambda=0;
    % beta=0;

    % while (iter<500 && norm(gra(:))>epsilon)  
       
    %     U_pre=U;
    %     V_pre=V;

    %     % backtracking   
    %     t=2;
    %     flag=false;
    %     for i=1:150
    %         t=rho*t;

    %         U=Uy-t*grad_U;
    %         V=Vy-t*grad_V;
    %         % fun_val_cur=fval(U,V,MA);
            
    %         % fprintf("pre function value is %3.4f, current function value is %3.4f, t is %3.4f\n",fval(U_pre,V_pre,MA), fun_val_cur, t);
    %         % fprintf("difference is %3.4f, tol is %3.4f\n",fval(U_pre,V_pre,MA)-fval(U,V,MA), alpha*t*norm(grad(:))^2);
    %         % gn=gradient(U,V,MA);
    %         % fprintf("norm of the gradient is %3.4f\n",norm(gn(:)));

    %         if fval(Uy,Vy,MA)-fval(U,V,MA)>=alpha*t*norm(gra(:))^2
    %             flag=true;
    %             fprintf("backtracking stepsize selection successful, step length is %3.4f\n",t);
    %             break;
    %         end
    %     end

    %     % update y
    %     % Uy=(1+beta)*U-beta*U_pre;
    %     % Vy=(1+beta)*V-beta*V_pre;
    %     beta=iter/(iter+3);
    %     Uy=U+beta*(U-U_pre);
    %     Vy=V+beta*(V-V_pre);


    %     % calculate the gradient
    %     gra=grad(Uy,Vy,MA);
    %     grad_U=gra(1:m,1:r);
    %     grad_V=gra(m+1:m+n,1:r);

    %     % update lambda
    %     lambda_pre=lambda;
    %     lambda=(1+sqrt(1+4*lambda^2))/2;
    %     fprintf("lambda is %3.8f\n",lambda);

    %     % update beta
    %     beta=(lambda_pre-1)/lambda;
    %     fprintf("beta is %3.8f\n",beta);
        
    %     fun_val=fval(Uy,Vy,MA);
    %     fprintf("iter_number=%d,norm_grad=%3.5f,fun_val=%3.5f\n",iter,norm(gra(:)),fun_val);

    %      iter=iter+1;
    % end  
    

    % initial point
    U=U0;
    V=V0;
    U_pre=U0;
    V_pre=V0;

    iter=0;

    [m,r]=size(U);
    [n,r]=size(V);
    
    fun_val=fval(U,V,MA); 
    fprintf("initial objective function value is %3.4f\n",fun_val);
    

    lambda=1;
   
    while (iter<1500) 
        iter=iter+1;
        % update lambda
        lambda_pre=lambda;
        lambda=(1+sqrt(1+4*lambda^2))/2;

        % update beta
        beta=(lambda_pre-1)/lambda;
        fprintf("beta is %3.8f\n",beta);

        Uy=U+beta*(U-U_pre);
        Vy=V+beta*(V-V_pre);

        % calculate the gradient
        gra=grad(Uy,Vy,MA);

        if norm(gra(:))<epsilon
            fprintf("============================== iter is %d, terminate ============================\n",iter);
            
        end    
        grad_U=gra(1:m,1:r);
        grad_V=gra(m+1:m+n,1:r);


        U_pre=U;
        V_pre=V;

        % line search
        t=2;
        flag=false;
        for i=1:150
            t=rho*t;

            U=Uy-t*grad_U;
            V=Vy-t*grad_V;
            % fun_val_cur=fval(U,V,MA);
            
            fprintf("pre function value is %3.10f, current function value is %3.10f, t is %3.10f\n",fval(Uy,Vy,MA), fval(U,V,MA), t);
            fprintf("difference is %3.10f, tol is %3.10f\n",fval(Uy,Vy,MA)-fval(U,V,MA), alpha*t*norm(gra(:))^2);
            % gn=gradient(U,V,MA);
            % fprintf("norm of the gradient is %3.4f\n",norm(gn(:)));

            if fval(Uy,Vy,MA)-fval(U,V,MA)>=alpha*t*norm(gra(:))^2
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.10f\n",t);
                break;
            end
        end        

        if ~flag
            fprintf("line search failed!\n");
        end    

        grad_U = grad(U,V,MA);


        if norm(grad_U(:))<epsilon
            fprintf("============================== iter is %d, terminate ============================\n",iter);
            break;
        end    
        fun_val=fval(U,V,MA);
        fprintf("iter_number=%d,norm_grad=%3.10f,norm_grad_U=%3.10f,fun_val=%3.5f\n",iter,norm(gra(:)),norm(grad_U(:)),fun_val);
    end    

    runhist.iter=iter;
    runhist.U=U;
    runhist.V=V;
    runhist.cput=cputime-timect;
end    