function runhist = Y2_subproblem (Aparas, B, cg_paras, paras)
    n=paras.n;
    r=paras.r;
    X = zeros(n-r,r);
    R=B;
    gamma=0;
    
    iter = 0;
    R_norm = norm(R(:));
   
    itmax=cg_paras.itmax;
    tol=cg_paras.tol;
    while (iter < itmax && R_norm>tol)
        %fprintf("X2_subproblem R norm is %3.4f\n", R_norm);
        % step 1
        iter = iter + 1;
       
        % step 2
        if iter == 1
            C = R;
        else
            beta = (norm(R(:))^2)/(norm(R_pre(:))^2);
            C=R+beta*C;
        end

        % Step 3
        alpha = norm (R(:))^2 / trace(C'*CadjointC(Aparas, paras, C,gamma));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*CadjointC(Aparas, paras, C,gamma);

        R_norm=norm(R(:));
    end

    %fprintf("X2_subproblem iter=%d,R norm is %3.12f\n",iter,R_norm);

    runhist.R_norm=R_norm;
    runhist.Y2=X;
    runhist.iter=iter;
    %fprintf("------------------ Y2_subproblem iter=%d ---------------------------\n",iter);
end
