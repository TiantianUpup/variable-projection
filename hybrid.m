function runhist = hybrid(U0,V0,MA,vp_paras,cg_paras,paras)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;

    % U_0 and V_0
    U = U0;
    V = V0;
    gra=grad(U,V,MA);

    fprintf("this is hybrid method, the initial objective value fval is %3.4f\n",fval(U,V,MA));

    iter = 0; 
    %itmax = paras.itmax;

    gamma=vp_paras.gamma;
    itmax=vp_paras.itmax;
    tol=vp_paras.tol;
    epsilon=1e-6;
    cg_paras.mu=vp_paras.mu;
    
    % equation (27)
    Ir=eye(r);
    J = zeros(r^2, r);
    J(1:r, 1:r) = Ir;  
    L = zeros(r^2, r^2-r);
    L(r+1:r^2, 1:r^2-r) = eye(r^2-r);
    
    % calculate the initial S
    R=proj_supp(U,V)*MA;
    I=eye((m+n)*r);
    Pk=10^-4*norm(R(:))*I;
    Ak=Pk;

    % calculate the permutation matrix and its blocks
    H = permutation(m,n,r);
    H1 = H(1:m*n,1:r^2);
    H21 = H(1:m*n,r^2+1:n*r);
    H22 = H(1:m*n,n*r+1:(n+m)*r-r^2);

    %&& norm(gra(:))>epsilon
    MA_temp=MA;
    while (iter < itmax && norm(gra(:))>tol)
        iter = iter + 1;
        % check whether (U\odot V)-\in K(U\odot V)
        UVkr=khatrirao(U,V);

        fprintf("the rank of the UVkr is %d\n",rank(UVkr));
        % UVkri=sym_inv(UVkr);
        % UViP=UVkri*UVkr*UVkri;
        % UV_res=UViP-UVkri;
        % fprintf("the residual is %3.50f\n", norm(UV_res(:))^2);

        % calculate tildeV_1 and tildeU_1 via the generalized QR decompositions of U and V

        [P, Sigma, S] = qr(U);
        S = S';
        [Q, Tau, T] = qr(V);
        T = T';
        Sigma_1 = Sigma(1:r, 1:r);
        tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
        Tau_1 = Tau(1:r, 1:r);
        tildeV_1 = Tau_1 * T;

        % fprintf("the norm of P is %3.10f\n",norm(P(:)));
        % fprintf("the norm of Q is %3.10f\n",norm(Q(:)));
        % fprintf("the norm of tildeU_1 is %3.10f\n",norm(tildeU_1(:)));
        % fprintf("the norm of tildeV_1 is %3.10f\n",norm(tildeV_1(:)));

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
        Aparas.S=Pk;
        
        % calculate B
        Z1=zeros(r,r);
        Z2=zeros(m-r,r);
        Z3=zeros(r,r);
        Z4=zeros(n-r,r);

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

        % Bparas.B1=B1;
        % Bparas.B2=B2;
        % Bparas.B3=B3;
        % paras.m = m;
        % paras.n = n;
        % paras.r = r;

        %% cg method
        runhist = cg_hy_subproblem(Aparas,cg_paras,paras,B0,gamma);

        %fprintf("%d-th variable projection method, X1_Y1_cg iter=%d, X2_cg iter=%d, Y2_cg iter=%d\n",iter, runhist_X1_Y1.iter,runhist_X2.iter,runhist_Y2.iter);

        x = runhist.x; % [X',Y']'
        X=reshape(x,m+n,r);

        % A=PX 
        Atemp=X(1:m,1:r);
        A=P*Atemp;
        % gauss_tol=1e-5;
        % if norm(A(:))<gauss_tol
        %     A=10*A;
        % end 

        fprintf("A norm is %3.10f\n",norm(A(:)));

        % B=QY   
        Btemp=X(m+1:m+n,1:r);
        B=Q*Btemp;
        fprintf("B norm is %3.10f\n",norm(B(:)));
        
        dir=[A;B];
        gra_gra_inner=-norm(gra(:))^2;
        gra_direction_inner=trace(gra'*dir);
        fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
    
        % if norm(dir(:))<epsilon
        %     fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
        %     break;
        % end    
    
        % if gra_direction_inner < gra_gra_inner 
        %     fprintf("-------------------------------- gauss-newton direction --------------------------------------------\n");
        % end    
            
        U_pre = U;
        V_pre = V;  
        fval_pre=fval(U_pre,V_pre,MA_temp);
        fprintf("fval(U_pre,V_pre,MA_temp) is %3.8f\n",fval_pre);
        
        flag = false;
        alpha=1e-10;
        t=1.25;
        beta=0.8;
        grad_norm=norm(gra(:))^2    ;
        for i=1:150
            t=beta*t;

            % U=U_pre+t*A;
            % V=V_pre+t*B;
            threshold=alpha*t*grad_norm;
            fval_cur=fval(U_pre+t*A,V_pre+t*B,MA_temp);
            fprintf("i is %d, fval(U_pre,V_pre,MA)-fval(U,V,MA) is %3.8f, alpha*t*norm(gra(:))^2 is %3.8f,t is %3.8f\n", i, fval_pre-fval_cur,threshold,t);
            
            % if fval_pre-fval_cur>0
            %     flag=true;
            %     fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
            %     break;
            % end    
            
            %if fval_pre-fval_cur>threshold
            if fval_pre-fval_cur>0
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
                break;
            end
        end

        if flag
            % update U
            U=U_pre+t*A;
            V=V_pre+t*B;
        else
            U=U_pre;
            V=V_pre;
            fprintf("line search failed!\n");
        end
        
        % update the gradient
        gra=grad(U,V,MA_temp);

        % update S and Ak
        R_pre=R;
        R=proj_supp(U,V)*MA_temp;
        RMA=R*MA_temp';
        Zk=(norm(R(:))/norm(R_pre(:)))*(gra+proj_diff_adj(U_pre,V_pre,RMA));
        zk=Zk(:);
        Sk=t*[A;B];
        sk=Sk(:);

        rela_fval=(fval(U_pre,V_pre,MA_temp)-fval(U,V,MA_temp))/fval(U_pre,V_pre,MA_temp);
        zs_inn=zk'*sk;
        fprintf("rela_fval is %3.10f, ZS_tra is %3.10f\n",rela_fval);
        fprintf("ZS_tra is %3.10f\n",zk'*sk);
        if zs_inn/norm(sk)^2>=1e-6
            fprintf("Approximation Hessian\n");
            AS=Ak*sk;
            Ak=Ak-(AS*AS')/(sk'*AS)+(zk*zk')/zs_inn;
            Sres=Ak*sk-zk;
            fprintf("rela_fval is %3.10f,Sres is %3.10f\n",zs_inn/norm(sk)^2,norm(Sres(:)));
           
            Pk=Ak;
        else
            Pk=1e-3*norm(R(:))*I;
            fprintf("simple choose, hybrid method, R_norm is %3.10f\n",norm(R(:)));
        end        

        % fprintf("the norm of the gauss-newton direction is %3.10f, epsilon*(norm(Sol(:))+epsilon) is %3.10f\n",t*norm(dir(:)),epsilon*(norm(Sol(:))+epsilon));
        % if t*norm(dir(:))<epsilon*(norm(Sol(:))+epsilon)
        %     break;
        % end  

        % if ~flag
        %     % gamma = 1;
        %     % U = U_pre + gamma * A;
        %     % V = V_pre + gamma * B;
        %     fprintf("line search failed!\n");
        % end

        fprintf("=========== the norm of the gradient is %3.8f ============\n",norm(gra(:)));

        fprintf("hybrid method %d-th iteration, fval is %3.20f =============================================================\n",iter,fval(U,V,MA_temp));
        
        fval_pre=fval(U_pre,V_pre,MA_temp);
        fval_cur=fval(U,V,MA_temp);
        if fval_pre-fval_cur<=1e-15*max(1,fval_cur)
            fprintf("the algorithm is terminate by the objective\n");
            break;
        end    
           fprintf("the norm of the gauss-newton direction is %3.10f\n",norm(dir(:)));
        if norm(dir(:))<epsilon
            fprintf("the algorithm is terminate by the gauss-newton direction condition\n");
            break;
        end
    end

    runhist.U = U;
    runhist.V = V;
    runhist.iter = iter;
    runhist.cput = cputime - timect;
end
