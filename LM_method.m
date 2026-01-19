function runhist = LM_method(U0,V0,MA,lm_paras,cg_paras,paras)
    fprintf("this is LM method\n");
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    % U_0 and V_0
    
    U=U0;
    V=V0;

    fprintf("the initial objective value fval is %3.4f\n",fval(U,V,MA));

    iter = 0;
    %itmax = paras.itmax;
    
    gamma=lm_paras.gamma;
    itmax=lm_paras.itmax;
    tol=lm_paras.tol; 
    epsilon=1e-6;

    % equation (27)
    J = zeros(r^2, r);
    J(1:r, 1:r) = eye(r);  
    L = zeros(r^2, r^2-r);
    L(r+1:r^2, 1:r^2-r) = eye(r^2-r);

    % calculate the permutation matrix and its blocks
    H = permutation(m,n,r);
    H1 = H(1:m*n,1:r^2);
    H21 = H(1:m*n,r^2+1:n*r);
    H22 = H(1:m*n,n*r+1:(n+m)*r-r^2);

    gra=grad(U,V,MA);

    %mu=lm_paras.mu;
    mu=1e-2;
    v=10;
    beta=lm_paras.beta;
    %nu=2;
    nu=beta;
    alpha=lm_paras.alpha;
    p=lm_paras.p;

    fprintf("============================== mu is %3.4f,gamma is %d,beta is %d,p is %d==================================\n",mu,alpha,beta,p);
    cg_paras.mu=mu;
    MA_temp=MA;

    % calculate tildeV_1 and tildeU_1 via the generalized QR decompositions of U and V
    [P, Sigma, S] = qr(U);
    S = S';
    [Q, Tau, T] = qr(V);
    T = T';
    Sigma_1 = Sigma(1:r, 1:r);
    tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
    Tau_1 = Tau(1:r, 1:r);
    tildeV_1 = Tau_1 * T;

    % calculate M via the generalized QR decomposition of G
    G = khatrirao(tildeU_1, tildeV_1);
    [W, Db, F] = qr(G);
    F = F';
    D = Db(1:r, 1:r);

    % calculate \overline{M(A)}
    kro = kron(P, Q);
    MA = kro' * MA_temp;
    
    % calculate M
    HWJ = H1 * W * J;
    Di=inv(D);% the dimension of D is r\times r
    M = F'*Di*HWJ'; 
    
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

    B0=[Z1;Z2;Z3;Z4];

    % the calculation of the initial mu
    while (iter < itmax && norm(gra(:))>tol)
        iter = iter + 1;

        %% cg method
        runhist = cg_subproblem(Aparas,cg_paras,paras,B0,gamma);

        %fprintf("%d-th variable projection method, X1_Y1_cg iter=%d, X2_cg iter=%d, Y2_cg iter=%d\n",iter, runhist_X1_Y1.iter,runhist_X2.iter,runhist_Y2.iter);

        X = runhist.X; % [X',Y']'
      
        % A=PX 
        Atemp=X(1:m,1:r);
        A=P*Atemp;
        % B=QY   
        Btemp=X(m+1:m+n,1:r);
        B=Q*Btemp;

        %% check whether [A;B] is the descent direction
        dir=[A;B];
        fprintf("the norm of the gauss-newton direction is %3.10f\n",norm(dir(:)));
        gra_gra_inner=-norm(gra(:))^2;
        gra_direction_inner=trace(gra'*dir);
        fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
    
        Sol=[U;V];
        % if norm(dir(:))<epsilon
        %     fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
        %     break;
        % end    
             
        % fprintf("the norm of the gauss-newton direction is %3.10f, epsilon*(norm(Sol(:))+epsilon) is %3.10f\n",norm(dir(:)),epsilon*(norm(Sol(:))+epsilon));
        % if norm(dir(:))<epsilon*(norm(Sol(:))+epsilon)
        %     break;
        % end 
        % fprintf("the norm of the gauss-newton direction is %3.10f\n",norm(dir(:)));

        fprintf("the norm of the gauss-newton direction is %3.10f\n",norm(dir(:)));
        if norm(dir(:))<epsilon
            fprintf("the algorithm is terminate by the gauss-newton direction condition\n");
            break;
        end

        % U_pre = U;
        % V_pre = V;  
        % fval_pre=fval(U_pre,V_pre,MA_temp);
        % fprintf("fval(U_pre,V_pre,MA_temp) is %3.8f\n",fval_pre);

        % calculate rho
        fval_pre=fval(U,V,MA_temp);
        fprintf("fval_pre is %3.4f\n",fval_pre);
        fval_cur=fval(U+A,V+B,MA_temp);
        fprintf("fval_cur is %3.4f\n",fval_cur);
    
        rho=fval_pre-fval_cur;
        %rho=(fval_pre-fval_cur)/Ldiff;
        % fprintf("Lgap is %3.10f\n",Lgap(U,V,A,B,MA_temp));
        if rho >0
            fprintf("rho is greater than 0 and rho is %3.4f\n",rho);
        else
            fprintf("rho is less than 0\n");
        end    

        if rho>0  
            % step acceptable, update U and V
            U=U+A;
            V=V+B;

            % % update mu and nu
            % mu=mu*max(1/3,1-(2*rho-1)^3);
            % nu=2;  

            % % method 1
            % mu=mu*max(1/alpha,1-(beta-1)*(2*rho-1)^p);
            % fprintf("rho is %3.10f, mu is %3.10f\n",rho,mu);
            % nu=beta;
            
            % method 2
            mu=mu/v;
            
            % update the gradient
            gra=grad(U,V,MA_temp);
            
            % update the A operator parameters
            % calculate tildeV_1 and tildeU_1 via the generalized QR decompositions of U and V
            [P, Sigma, S] = qr(U);
            S = S';
            [Q, Tau, T] = qr(V);
            T = T';
            Sigma_1 = Sigma(1:r, 1:r);
            tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
            Tau_1 = Tau(1:r, 1:r);
            tildeV_1 = Tau_1 * T;

            % calculate M via the generalized QR decomposition of G
            G = khatrirao(tildeU_1, tildeV_1);
            [W, Db, F] = qr(G);
            F = F';
            D = Db(1:r, 1:r);

            % calculate \overline{M(A)}
            kro = kron(P, Q);
            MA = kro' * MA_temp;
            
            % calculate M
            HWJ = H1 * W * J;
            Di=inv(D);% the dimension of D is r\times r
            M = F'*Di*HWJ'; 
            
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

            B0=[Z1;Z2;Z3;Z4];

            fprintf("iter=%d, fval is %3.8f\n",iter,fval_cur);
        else
            % % update mu and nu (method 1)
            % mu=mu*nu;
            % nu=2*nu;

            % method 2
            mu=mu*v;

            fprintf("iter=%d, fval is %3.8f\n",iter,fval_pre);
        end    


        % update the cg parameter
        cg_paras.mu=mu;
        fprintf("============================================= mu is %3.10f ============================\n",mu);
       
        fprintf("=========== the norm of the gradient is %3.8f ============\n",norm(gra(:)));

        fprintf("LM method %d-th iteration, fval is %3.20f =============================================================\n",iter,fval(U,V,MA_temp));        
    end

    runhist.U = U;
    runhist.V = V;
    runhist.iter = iter;
    runhist.cput = cputime - timect;
end
