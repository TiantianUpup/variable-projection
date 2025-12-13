function runhist = variable_projection(MA, rho, paras,X)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    % U_0 and V_0
    U = rand(m,r);
    V = rand(n,r);

    
    fprintf("the initial objective value fval is %3.4f\n",fval(U,V,MA));

    iter = 1;
    %itmax = paras.itmax;
    itmax=10;

    % equation (27)
    J = zeros(r^2, r);
    J(1:r, 1:r) = eye(r);  
    L = zeros(r^2, r^2-r);
    L(r+1:r^2, 1:r^2-r) = eye(r^2-r);

    MA_temp=MA;
    while (iter <= itmax)
        % check whether (U\odot V)-\in K(U\odot V)
        % UVkr=khatrirao(U,V);
        % fprintf("the rank of the UVkr is %d\n",rank(UVkr));
        % [a,b]=size(UVkr);
        % UVkrr=rank(UVkr);
        % [A,B,C]=qr(UVkr);
        % A11=B(1:UVkrr,1:UVkrr);
        % Z=zeros(b,a);
        % Z(1:UVkrr,1:UVkrr)=inv(A11);
        % UVi=C*Z*A';
        % UViP=UVi*UVkr*UVi;
        % UV_res=UViP-UVi;
        % fprintf("the residual is %3.8f\n", norm(UV_res(:))^2);

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
        MA = kro' * MA;

        % calculate the permutation matrix and its blocks
        H = permutation(m,n,r);
        H1 = H(1:m*n,1:r^2);
        H21 = H(1:m*n,r^2+1:n*r);
        H22 = H(1:m*n,n*r+1:(n+m)*r-r^2);
        
        % calculate M via the generalized QR decomposition of G
        G = khatrirao(Sigma_1 * S, Tau_1 * T);
        [W, Db, F] = qr(G);
        F = F';
        D = Db(1:r, 1:r);

        HWJ = H1 * W * J;
        % inv(D) = D\I
        M = F' * inv(D) * HWJ';
        %[l,~]=size(D);
        %fprintf("the dimension of D is %d, and the rank of D is %d\n",l,rank(D));
        

        % calculate B1
        WL = W * L;
        B1 = WL' * H1' * MA;
        % calculate B2
        B2 = H22' * MA;
       
        % calculate B3
        B3 = H21' * MA;

        Aparas.tildeU_1 = tildeU_1;
        Aparas.tildeV_1 = tildeV_1;
        Aparas.M = M;
        Aparas.MA = MA;
        Aparas.WL = WL;

        % paras.m = m;
        % paras.n = n;
        % paras.r = r;

        %% cg method
        runhist_X1_Y1 = X1_Y1_subproblem(Aparas, B1, paras);
        runhist_X2 = X2_subproblem(Aparas, B2, paras);
        runhist_Y2 = Y2_subproblem(Aparas, B3, paras);

        %fprintf("%d-th variable projection method, X1_Y1_cg iter=%d, X2_cg iter=%d, Y2_cg iter=%d\n",iter, runhist_X1_Y1.iter,runhist_X2.iter,runhist_Y2.iter);

        Temp = runhist_X1_Y1.X; % Temp=[X_1',Y_1']'
        % fprintf("Temp is\n");
        % disp(Temp);
        X2 = runhist_X2.X2;
        Y2 = runhist_Y2.Y2;

        A = zeros(m, r);
        B = zeros(n, r);
        A(1:r,1:r) = Temp(1:r, 1:r);
        A(r+1:m,1:r) = X2;
        B(1:r,1:r) = Temp(r + 1:2 * r, 1:r);
        B(r+1:n,1:r) = Y2;

        % check whether (A,B) is the descent direction
        GN_direction=[A;B];
        gra=gradient(U,V,MA_temp);
        gra_direction_inner=trace(gra'*GN_direction);
        fprintf("========== gra_direction_inner is %3.4f\n",gra_direction_inner);

        % fprintf("A is \n");
        % disp(A);
        % fprintf("B is \n");
        % disp(B);

        % U_pre=U;
        % V_pre=V;
        
        % fprintf("the residual of U and U_pre is %3.4f\n",norm(res(:))^2);

        % fprintf("================== U is =========================\n");
        % disp(U);
        % fprintf("================== V is =========================\n");
        % disp(V);

        %% line search procedure
        gamma = 1;
        U_pre = U;
        V_pre = V;
        flag = false;
 
        k=0;
        for i = 1:100
            k=k+1; 
            gamma = gamma * rho;
    
            U = U_pre + gamma * A;
            V = V_pre + gamma * B;
  
            %fprintf("fval(U, V,MA_temp) is %3.8f,fval(U_pre, V_pre,MA_temp) %3.8f, gamma is %3.8f\n",fval(U, V,MA_temp),fval(U_pre, V_pre,MA_temp),gamma);

            if (fval(U, V,MA_temp) <= (1 - gamma) * fval(U_pre, V_pre,MA_temp))
                flag = true;
                fprintf('k=%d, gamma is %3.8f, line search successful!\n',k,gamma);
                break;
            end

        end

        if ~flag
            gamma = 0.001;
            U = U_pre + gamma * A;
            V = V_pre + gamma * B;
            fprintf('line search failed!\n');
        end

        kr=khatrirao(U,V);
        krmp=pinv(kr); % the Moore-Penrose inverse of U\odot V
        W=MA_temp'*krmp';
        Xhat=generate_cp_tensor(U,V,W);
        Xres=Xhat-X;

        fprintf("variable projection %d-th iteration, fval is %3.4f, residual is %3.4f\n",iter,fval(U,V,MA_temp),norm(Xres(:))^2);
        iter = iter + 1;
    end

    runhist.U = U;
    runhist.V = V;
    runhist.iter = iter;
    runhist.vp_cput = cputime - timect;
end
