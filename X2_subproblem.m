function runhist = X2_subproblem (Aparas, B, paras)
    % parameters of the operator A
    tildeV_1 = Aparas.tildeV_1;
    M = Aparas.M;
    MA = Aparas.MA;
    MMA = M*MA;

    % convergence condition
    %itmax=paras.itmax;
    itmax=150;
    %tol=paras.tol;
    tol=1e-4;

    % parameter dimension
    m=paras.m;
    r=paras.r;

    % initialization
    X = zeros(m-r,r);
    % calculation of R
    R = zeros(m-r,r);
    I = eye(m-r);
    BM = B*MMA';
   
    for i=1:r
        tildeV_1_i=tildeV_1(:,i);
        R(:,i)=kron(I,tildeV_1_i')*BM(:,i);
    end 

    iter = 0;
    R_norm = norm(R(:));
   
    
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
        alpha = norm (R(:))^2 / trace(C'*CadjointC(tildeV_1, MMA, paras,C));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*CadjointC(tildeV_1, MMA, paras,C);

        R_norm=norm(R(:));
    end

    %fprintf("X2_subproblem iter=%d,R norm is %3.12f\n",iter,R_norm);

    runhist.R_norm=R_norm;
    runhist.X2=X;
    runhist.iter=iter;
end
