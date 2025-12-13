function runhist = X1_Y1_subproblem (Aparas, B, paras)
    % parameters of the operator A
    tildeV_1 = Aparas.tildeV_1;
    tildeU_1 = Aparas.tildeU_1;
    M = Aparas.M;
    MA = Aparas.MA;
    WL=Aparas.WL;
    MMA = M*MA;

    % convergence condition
    %itmax=paras.itmax;
    itmax=150;
    %tol=paras.tol;
    tol=1e-4;

    % initialization
    r=paras.r;
    X = zeros(2*r,r);
    % calculation of R
    R = zeros(2*r,r);
    R1=zeros(r,r);
    R2=zeros(r,r);
    I = eye(r);

    % calculate R
    WBS=WL*B*MMA'; 
    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        tildeV_1_i=tildeV_1(:,i);
        R1(:,i)=kron(I,tildeV_1_i')*WBS(:,i);
        R2(:,i)=kron(tildeU_1_i',I)*WBS(:,i);
    end                                                          

    R(1:r,1:r)=R1;
    R(r+1:2*r,1:r)=R2;

    iter = 0;
    R_norm = norm(R(:));
    %fprintf("X1_Y1_subproblem R norm is %3.4f\n", R_norm);

    while (iter  < itmax && R_norm>tol)
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
        alpha = norm (R(:))^2 / trace(C'*AadjointA(tildeU_1,tildeV_1,MMA,WL,paras,C));

        % Step 4
        X = X + alpha * C;

        % Step 5
        R_pre=R;
        R = R - alpha*AadjointA(tildeU_1,tildeV_1,MMA,WL,paras,C);

        R_norm=norm(R(:));
    end

    %fprintf("X1_Y1_subproblem iter=%d\n",iter);

    runhist.R_norm=R_norm;
    runhist.X=X;
    runhist.iter=iter;
end    