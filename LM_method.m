function runhist = LM_method(MA, rho,gamma,paras)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    % U_0 and V_0
    U = rand(m,r);
    V = rand(n,r);

    U0=U;
    V0=V;

    fprintf("the initial objective value fval is %3.4f\n",fval(U,V,MA));

    iter = 1;
    %itmax = paras.itmax;
    itmax=150;
    epsilon=1e-4;
    % equation (27)
    J = zeros(r^2, r);
    J(1:r, 1:r) = eye(r);  
    L = zeros(r^2, r^2-r);
    L(r+1:r^2, 1:r^2-r) = eye(r^2-r);

    gra=grad(U,V,MA);
    lambda=1e-3;
    beta_1=0.3;
    beta_2=2;

    %&& norm(gra(:))>epsilon
    MA_temp=MA;
    while (iter <= itmax)
        % calculate tildeV_1 and tildeU_1 via the generalized QR decompositions of U and V
        [P, Sigma, S] = qr(U);
        S = S';
        [Q, Tau, T] = qr(V);
        T = T';
        Sigma_1 = Sigma(1:r, 1:r);
        tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
        Tau_1 = Tau(1:r, 1:r);
        tildeV_1 = Tau_1 * T;

        % calculate \overline{M(A)}
        kro = kron(P, Q);
        MA = kro' * MA_temp;

        % calculate the permutation matrix and its blocks
        H = permutation(m,n,r);
        H1 = H(1:m*n,1:r^2);
        H21 = H(1:m*n,r^2+1:n*r);
        H22 = H(1:m*n,n*r+1:(n+m)*r-r^2);
        
        % calculate M via the generalized QR decomposition of G
        G = khatrirao(Sigma_1 * S, Tau_1 * T);
        [W, Db, F] = qr(G);
        %F = F';
        D = Db(1:r, 1:r);

        HWJ = H1 * W * J;
        % inv(D) = D\I
        Di=inv(D);
        M = F*Di*HWJ'; 
        %[l,~]=size(D);
        %fprintf("the dimension of D is %d, and the rank of D is %d\n",l,rank(D));
        
        WL = W * L;
        Aparas.A=WL;
        MMA=M*MA;
        Aparas.C=MMA;
        Aparas.E=Di'*F;
        Aparas.tildeU_1=tildeU_1;
        Aparas.tildeV_1=tildeV_1;
        B1=WL'*H1'*MA;
        B2=H21'*MA;
        B3=H22'*MA;
        Aparas.G=WL*B1;
        Aparas.B2=B2;
        Aparas.B3=B3;
        
        % calculate B
        Z1=zeros(r,r);
        Z2=zeros(m-r,r);
        Z3=zeros(r,r);
        Z4=zeros(n-r,r);

        Ir=eye(r);
        Im=eye(m-r);
        In=eye(n-r);

        Temp1=WL*B1*MMA';
        Temp2=B3*MMA';
        Temp3=B2*MMA';
        for i=1:r
            tildeU_1_i=tildeU_1(:,i);
            tildeV_1_i=tildeV_1(:,i);
            Z1(:,i)=kron(Ir,tildeV_1_i')*Temp1(:,i);
            Z2(:,i)=kron(Im,tildeV_1_i')*Temp2(:,i);
            Z3(:,i)=kron(tildeU_1_i',Ir)*Temp1(:,i);
            Z4(:,i)=kron(tildeU_1_i',In)*Temp3(:,i);
        end   

        B=[Z1;Z2;Z3;Z4];

        %% cg method
        runhist = cg_subproblem(Aparas,paras,B,gamma,lambda);

        
        X = runhist.X; % [X',Y']'
        % A=PX 
        Atemp=X(1:m,1:r);
        A=P*Atemp;
        % B=QY   
        Btemp=X(m+1:m+n,1:r);
        B=Q*Btemp;


        %% check whether [A;B] is the descent direction
        dir=[A;B];
        gra=grad(U,V,MA_temp);
        gra_direction_inner=trace(gra'*dir);
        fprintf("gra_direction_inner is %3.4f\n",gra_direction_inner);
    
        % f=fval(U+A,V+B,MA_temp);
        % f_pre=fval(U,V,MA_temp);
        % fprintf("fval(U+A,V+B,MA_temp) is %3.4f, fval(U,V,MA_temp) is %3.4f\n", f,f_pre);
        fprintf("fval(U+A,V+B,MA_temp) is %3.4f, fval(U,V,MA_temp) is %3.4f\n", fval(U+A,V+B,MA_temp),fval(U,V,MA_temp));
        if fval(U+A,V+B,MA_temp)<fval(U,V,MA_temp)
            fprintf("(A,B) is a decrease direction\n");
            U=U+A;
            V=V+B;
            % decrease lambda
            lambda=beta_1*lambda;
        else
            fprintf("(A,B) is not a decrease direction\n");
            % increase lambda
            lambda=beta_2*lambda;
        end         

        % update the gradient
        gra=grad(U,V,MA_temp);
        fprintf("=========== the norm of the gradient is %3.8f ============\n",norm(gra(:)));

        fprintf("variable projection %d-th iteration, fval is %3.8f\n",iter,fval(U,V,MA_temp));
        iter = iter + 1;
    end

    runhist.U = U;
    runhist.V = V;
    runhist.iter = iter;
    runhist.vp_cput = cputime - timect;
end
