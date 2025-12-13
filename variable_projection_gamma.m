function runhist = variable_projection_gamma(MA, rho,gamma,paras)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    % U_0 and V_0
    U = rand(m,r);
    V = rand(n,r);

    U0=U;
    V0=V;

    fprintf("variable_projection_gamma method, the initial objective value fval is %3.4f\n",fval(U,V,MA));

    iter = 1; 
    %itmax = paras.itmax;
    itmax=15;
    epsilon=1e-4;
    lambda=0;
    % equation (27)
    J = zeros(r^2, r);
    J(1:r, 1:r) = eye(r);  
    L = zeros(r^2, r^2-r);
    L(r+1:r^2, 1:r^2-r) = eye(r^2-r);

    gra=grad(U,V,MA);

    %&& norm(gra(:))>epsilon
    MA_temp=MA;
    while (iter <= itmax && norm(gra(:))>epsilon)
        % % check whether (U\odot V)-\in K(U\odot V)
        % UVkr=khatrirao(U,V);
        % fprintf("the rank of the UVkr is %d\n",rank(UVkr));
        % [a,b]=size(UVkr);
        % UVkrr=rank(UVkr);
        % [A,B,C]=qr(UVkr);
        % A11=B(1:UVkrr,1:UVkrr);
        % Z=zeros(b,a);
        % Z(1:UVkrr,1:UVkrr)=inv(A11);
        % UVi=C'*Z*A';
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

        % Bparas.B1=B1;
        % Bparas.B2=B2;
        % Bparas.B3=B3;
        % paras.m = m;
        % paras.n = n;
        % paras.r = r;

        %% cg method
        runhist = cg_subproblem(Aparas,paras,B,gamma,lambda);

        %fprintf("%d-th variable projection method, X1_Y1_cg iter=%d, X2_cg iter=%d, Y2_cg iter=%d\n",iter, runhist_X1_Y1.iter,runhist_X2.iter,runhist_Y2.iter);

        X = runhist.X; % [X',Y']'
       
        % A = zeros(m, r);
        % B = zeros(n, r);
        % A(1:r,1:r) = X(1:r,1:r);
        % A(r+1:m,1:r) = X(r+1:m,1:r);
        % B(1:r,1:r) = X(m+1:m+r,1:r);
        % B(r+1:n,1:r) = X(m+r+1:m+n,1:r);

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
        %if gra_direction_inner>0
            value_temp=proj_supp(U,V)*(khatrirao(U,B)+khatrirao(A,V))*sym_inv(khatrirao(U,V));
            value=value_temp+value_temp';

            fprintf("when inner is greater than 0, Jp is %3.8f\n",norm(value(:)));
        %end

        % if gra_direction_inner<0
        %     fprintf("gauss-newton direction, the norm of the gauss-newton direction is %3.8f\n",norm(dir(:)));
        % else
        %     A=-gra(1:m,1:r);
        %     B=-gra(m+1:m+n,1:r);
        %     dir=[A;B];

        %     fprintf("negative gradient direction, the norm of the negative gradient direction is %3.8f\n",norm(dir(:)));
        %     % gra=grad(U,V,MA_temp);
        %     % gra_direction_inner=trace(gra'*dir);
        %     %fprintf("gra_direction_inner is %3.4f\n",gra_direction_inner);
        % end        

        U_pre = U;
        V_pre = V;  
        fprintf("fval(U_pre,V_pre,MA_temp) is %3.8f\n",fval(U_pre, V_pre,MA_temp));

        %% line search procedure
        %gamma = 1;
        % U_pre = U;
        % V_pre = V;
        flag = false;
   
        % fprintf("fval(U_pre, V_pre,MA_temp) is %3.8f\n",fval(U_pre, V_pre,MA_temp));
        % for i = 1:100
        %     gamma = gamma * rho;
    
        %     U = U_pre + gamma * A;
        %     V = V_pre + gamma * B;
  
        %     fprintf("fval(U, V,MA_temp) is %3.8f, (1-gamma)*fval(U_pre, V_pre,MA_temp) %3.8f, gamma is %3.8f\n",fval(U, V,MA_temp), (1-gamma)*fval(U_pre, V_pre,MA_temp),gamma);

        %     if (fval(U, V,MA_temp) <= (1 - gamma) * fval(U_pre, V_pre,MA_temp))
        %         flag = true;
        %         fprintf('line search successful!\n');
        %         break;
        %     end

        % end

        % if ~flag
        %     gamma = 1;
        %     U = U_pre + gamma * A;
        %     V = V_pre + gamma * B;
        %     fprintf('line search failed!\n');
        % end

        alpha=1e-3;
        t=1;
        beta=0.8;
        for i=1:150
            t=beta*t;

            U=U_pre+t*A;
            V=V_pre+t*B;
            
            fprintf("i is %d, fval(U_pre,V_pre,MA)-fval(U,V,MA) is %3.8f, alpha*t*norm(gra(:))^2 is %3.8f,t is %3.8f\n", i, fval(U_pre,V_pre,MA_temp)-fval(U,V,MA_temp), alpha*t*norm(gra(:))^2,t);
            if fval(U_pre,V_pre,MA_temp)-fval(U,V,MA_temp)>alpha*t*norm(gra(:))^2
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
                break;
            end
        end

        if ~flag
            % gamma = 1;
            % U = U_pre + gamma * A;
            % V = V_pre + gamma * B;
            fprintf("line search failed!\n");
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
