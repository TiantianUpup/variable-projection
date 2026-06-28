function runhist = tcg_subproblem(Aparas,cg_paras,paras,B,delta)
    timect = cputime;

    % convergence condition
    itmax=cg_paras.itmax;
    tol=cg_paras.tol;
    
    % initialization
    r=paras.r;
    m=paras.m;
    n=paras.n;
    
    Z=zeros(m+n,r);
    R=-B;
    D=B;

    iter=0; 

    while (iter < itmax)
        if iter == 0 && norm(R(:))<tol
            P=zeros(m+n,r);
            break;
        end 

        iter = iter + 1;
        
        Bd = Aope_simple(Aparas, paras, D);
        dTBd=sum(D(:).*Bd(:));

        if dTBd <= 0
            a = sum(D(:).^2);
            b = sum(Z(:).*D(:));
            c = sum(Z(:).^2)-delta^2;
            tau =(-b+sqrt(b^2-a*c))/a;

            fprintf("-------------- iter=%d, dTBd<=0, tau=%d --------------\n",iter,tau);

            P = Z + tau*D;
            break;
        end 
        
        % update alpha
        alpha = sum(R(:).^2)/dTBd;

        %update z
        Z = Z+alpha*D;

        if norm(Z(:)) >= delta
            a = sum(D(:).^2);
            b = sum(Z(:).*D(:));
            c = sum(Z(:).^2)-delta^2;
            tau =(-b+sqrt(b^2-a*c))/a;

             fprintf("-------------- iter=%d, Znorm>=delta, tau=%d --------------\n",iter,tau);

            P = Z + tau*D;
            break;
        end   
        
        % update 
        R_pre = R;
        R = R + alpha*Aope_simple(Aparas, paras, D);
        R_norm = norm(R(:));
        if R_norm < tol
            P=Z;
            break;
        end    

        % update beta
        beta = sum(R(:).^2)/sum(R_pre(:).^2);
        
        % update D
        D=-R + beta*D;
    end

    cput = cputime - timect;
    fprintf("cg method iter=%d, R_norm=%3.7f, time cost is %3.4f\n",iter,norm(R_norm(:)),cput);
    runhist.R_norm = R_norm;
    runhist.X = P;
    runhist.iter = iter;
end
