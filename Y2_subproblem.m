function runhist = Y2_subproblem (Aparas, B, paras)
    % parameters of the operator A
    tildeU_1 = Aparas.tildeU_1;
    M = Aparas.M;
    MA = Aparas.MA;
    MMA = M*MA;

    % convergence condition
    %itmax=paras.itmax;
    itmax=150;
    %tol=paras.tol;
    tol=1e-4;

    % parameter dimension
    n=paras.n;
    r=paras.r;

    % initialization
    Y2 = zeros(n-r,r);
    % calculation of R
    R = zeros(n-r,r);
    I = eye(n-r);
    MMAt=MMA';
    
    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        R(:,i)=kron(tildeU_1_i',I)*B*MMAt(:,i);
    end 

    iter = 0;
    R_norm = norm(R(:));

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
        alpha = norm(R(:))^2 / trace(C'*BadjointB(tildeU_1, MMA, paras,C));

        % Step 4
        Y2 = Y2 + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*BadjointB(tildeU_1, MMA, paras,C);

        R_norm=norm(R(:));
    end

    %fprintf("Y2_subproblem iter=%d, R_norm is %3.12f\n",iter, R_norm);
    runhist.R_norm=R_norm;
    runhist.Y2=Y2;
    runhist.iter=iter;
end    