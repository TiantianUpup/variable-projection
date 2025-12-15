function [U,V] = gradient_descent(U0,V0,MA,alpha,beta,epsilon,X)
    % initial point
    U=U0;
    V=V0;
    gra = grad(U,V,MA);
    [m,r]=size(U);
    [n,r]=size(V);
    grad_U=gra(1:m,1:r);
    grad_V=gra(m+1:m+n,1:r);
  
    fun_val=fval(U,V,MA); 
    fprintf("initial objective function value is %3.4f\n",fun_val);
    iter=0;


    while (iter<5000 && norm(gra(:))>epsilon)
        kr=khatrirao(U,V);
        fprintf("the rank of kr is %d\n", rank(kr));
        
        iter=iter+1;
        
        % backtracking
        t=2;
        % while(fun_val-fval(U-t*grad_U,V-t*grad_V,MA)<alpha*t*norm(grad)^2)
        %     t=beta*t;
        % end    
        U_pre=U;
        V_pre=V;
        flag=false;
        for i=1:150
            t=beta*t;

            U=U_pre-t*grad_U;
            V=V_pre-t*grad_V;
            % fun_val_cur=fval(U,V,MA);
            
            % fprintf("pre function value is %3.4f, current function value is %3.4f, t is %3.4f\n",fval(U_pre,V_pre,MA), fun_val_cur, t);
            % fprintf("difference is %3.4f, tol is %3.4f\n",fval(U_pre,V_pre,MA)-fval(U,V,MA), alpha*t*norm(grad(:))^2);
            % gn=gradient(U,V,MA);
            % fprintf("norm of the gradient is %3.4f\n",norm(gn(:)));

            if fval(U_pre,V_pre,MA)-fval(U,V,MA)>=alpha*t*norm(gra(:))^2
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.4f\n",t);
                break;
            end


            % fprintf("pre function value is %3.4f, current function value is %3.4f, t is %3.4f\n",fun_val, fun_val_cur, t);
            % fprintf("(1-t) * fval(U_pre, V_pre,MA) is %3.4f\n", (1-t) * fval(U_pre, V_pre,MA));
            % if (fval(U, V,MA) <= (1-t) * fval(U_pre, V_pre,MA))
            %     flag=true;
            %     fprintf("backtracking stepsize selection successful, step length is %3.4f\n",t);
            %     break;
            % end    
        end

        % if ~flag
        %     t=1;
        %     U = U-t*grad_U;
        %     V = V-t*grad_V;
        %     fprintf('line search failed!\n');
        % end    
      
        % if norm(gra(:))<epsilon
        %     t=1;
        %     U = U-t*grad_U;
        %     V = V-t*grad_V;

        %     fprintf("escape from the saddle point!\n");
        % end    

        fun_val=fval(U,V,MA);

        % update the gradient
        gra=grad(U,V,MA);
        grad_U=gra(1:m,1:r);
        grad_V=gra(m+1:m+n,1:r);


        kr=khatrirao(U,V);
        krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
        What=MA'*krmp';
        Xhat=generate_cp_tensor(U,V,What);
        res=X-Xhat;

        fprintf("iter_number=%d,norm_grad=%3.5f,fun_val=%3.5f,res=%3.5f\n",iter,norm(gra(:)),fun_val,norm(res(:)));
    end    
end    