function runhist = X2_subproblem (Aparas, B, paras)
    gamma=0;

    % parameter dimension
    m=paras.m;
    r=paras.r;

    % initialization
    X2 = zeros(m-r,r);
    R=B;

    iter = 0;
    R_norm = norm(R(:));

    itmax=500;
    tol=1e-6;
    while (iter < itmax && R_norm>tol)
        %fprintf("Y2_subproblem R norm is %3.4f\n", R_norm);
        % step 1
        iter = iter + 1;

        % step 2
        if iter == 1
            C = R;
        else
            beta = norm(R(:))^2/norm(R_pre(:))^2;
            C=R+beta*C;
        end

        % Step 3
        alpha = norm(R(:))^2 / trace(C'*BadjointB(Aparas, paras, C,gamma));

        % Step 4
        X2 = X2 + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*BadjointB(Aparas, paras, C,gamma);

        R_norm=norm(R(:));
    end

    %fprintf("Y2_subproblem iter=%d, R_norm is %3.12f\n",iter, R_norm);
    runhist.R_norm=R_norm;
    runhist.X2=X2;
    runhist.iter=iter;
    %fprintf("------------------ X2_subproblem iter=%d ---------------------------\n",iter);
end    