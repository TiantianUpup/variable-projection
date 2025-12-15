function [U,V] = gauss_newton_original(MA, paras)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    % U_0 and V_0
    U = rand(m,r);
    V = rand(n,r);

    fprintf("variable_projection_gamma method, the initial objective value fval is %3.4f\n",fval(U,V,MA));

    % calculate B
    %B=-grad(U,V,MA);
    gra=grad(U,V,MA);
    iter=0;

    itmax=150;
    epsilon=1e-4;
    while (iter < itmax && norm(gra(:))>epsilon)
        iter=iter+1;
        Lparas.U=U;
        Lparas.V=V;
        Lparas.MA=MA;

        runhist = gn_cg_original (Lparas, paras, -gra);
        X=runhist.X;

        A=X(1:m,1:r);
        B=X(m+1:m+n,1:r);
        gra_direction_inner=trace(gra'*X);
        fprintf("gra_direction_inner is %3.4f\n",gra_direction_inner);
    
        U_pre = U;
        V_pre = V;  
        fprintf("fval(U_pre,V_pre,MA) is %3.8f\n",fval(U_pre, V_pre,MA));

        alpha=1e-3;
        t=1;
        beta=0.8;
        for i=1:150
            t=beta*t;

            U=U_pre+t*A;
            V=V_pre+t*B;
            fprintf("i is %d, fval(U,V,MA) is %3.8f\n", i, fval(U,V,MA));
            fprintf("i is %d, fval(U_pre,V_pre,MA)-fval(U,V,MA) is %3.8f, alpha*t*norm(gra(:))^2 is %3.8f,t is %3.8f\n", i, fval(U_pre,V_pre,MA)-fval(U,V,MA), alpha*t*norm(gra(:))^2,t);
            if fval(U_pre,V_pre,MA)-fval(U,V,MA)>alpha*t*norm(gra(:))^2
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
                break;
            end
        end

        % update the gradient
        gra=grad(U,V,MA);
        fprintf("iter=%d, gradient norm is %3.8f\n", iter, norm(gra(:)));
    end

    
end    