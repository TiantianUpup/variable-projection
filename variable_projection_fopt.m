function runhist = variable_projection_fopt(U0,V0,MA,vp_paras,cg_paras,paras)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;

    % U_0 and V_0
    U = U0;
    V = V0;

    [P, Sigma, S] = qr(U);
    [Q, Tau, T] = qr(V);
    S = S';
    T = T';

    P_tru=P(:,1:r);
    Sigma_tru=Sigma(1:r,:);

    Q_tru=Q(:,1:r);
    Tau_tru=Tau(1:r,:);
    X_hat=Sigma_tru*S;
    Y_hat=Tau_tru*T;
    XYkr=khatrirao(X_hat,Y_hat);

    [P_hat,Sigma_hat,Q_hat]=svd(XYkr,"econ");

    PQkr=kron(P_tru,Q_tru);
    PQP=PQkr*P_hat;
    Proj=PQP*PQP';

    % update the fval 
    R=MA-Proj*MA;
    fval_pre_opt=0.5*norm(R(:))^2;
    fval_init_opt=fval_pre_opt;

    % update the gradient
    UVi=Q_hat'*inv(Sigma_hat)*P_hat'*PQkr';
    UVi_qr=sym_inv(khatrirao(U,V));
    UVi_res=UVi-UVi_qr;
    fprintf("--------------- the residual of the UVi is %3.20f -----------\n",norm(UVi_res));

    PUVi=PQkr*P_hat*inv(Sigma_hat)*Q_hat;
    gra_opt_opt=grad_opt_opt_opt(R,UVi,PUVi,U,V,MA);

    % method 1
    % gra=grad(U,V,MA);
    UVkr=khatrirao(U,V);
    Pro=proj_supp(UVkr);
    PMA=Pro*MA;
    PMA_res=PMA-R;
    fprintf("--------------------------------- the residual of the PMA_res is %3.10f-----------------\n",norm(PMA_res(:)));
    gra=grad_opt(PMA,UVkr,U,V,MA);
    grad_res=gra_opt-gra;
    fprintf("the residual of grad is %3.20f\n",norm(grad_res(:)));

    %fval_pre=fval(U,V,MA);
    fval_pre=fval_opt(PMA);
    fval_init=fval_pre;

    fval_res=fval_pre-fval_pre_opt;
    grad_res2=gra_opt_opt-gra;


    % calculate the Z
    RMA=R*MA';
    %Z1=RMA*UVi'+RMA*PUVi;
    Z1=(PMA*MA'+PMA*PMA')*UVi_qr';
    Z2=(PMA*MA'+PMA*PMA')*UVi';
    Z_res=Z1-Z2;

    fprintf("--------------------------- the residual of Z is %3.20f --------------------------\n",norm(Z_res(:)));


    fprintf("-------------- the residual of f is %3.20f, the residual of the gradient is %3.20f--------------------------------\n",fval_res,norm(grad_res2(:)));

    fprintf("variable_projection_fopt method, the initial objective value fval is %3.4f\n",fval(U,V,MA));

    iter = 0; 
    %itmax = paras.itmax;

    gamma=vp_paras.gamma;
    itmax=vp_paras.itmax;
    tol=vp_paras.tol; 
    epsilon=1e-6;
    cg_paras.mu=vp_paras.mu;
    
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

    %&& norm(gra(:))>epsilon
    MA_temp=MA;
    while (iter < itmax && norm(gra(:))>tol)
        iter = iter + 1;
        % check whether (U\odot V)-\in K(U\odot V)
        %UVkr=khatrirao(U,V);

        %fprintf("the rank of the UVkr is %d\n",rank(UVkr));
        % UVkri=sym_inv(UVkr);
        % UViP=UVkri*UVkr*UVkri;
        % UV_res=UViP-UVkri;
        % fprintf("the residual is %3.50f\n", norm(UV_res(:))^2);

        % calculate tildeV_1 and tildeU_1 via the generalized QR decompositions of U and V

        
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
        
        % calculate B
        Z1=zeros(r,r);
        Z2=zeros(m-r,r);
        Z3=zeros(r,r);
        Z4=zeros(n-r,r);

        Temp1=WL*B1*MMA';
        Temp2=B3*MMA';
        Temp3=B2*MMA';

        for i=1:r
            tildeU_1_i=tildeU_1(:,i);
            tildeV_1_i=tildeV_1(:,i);

            Temp1_mat=reshape(Temp1(:,i),[r,r]);
            Temp2_mat=reshape(Temp2(:,i),[r,m-r]);

            Z1(:,i)=Temp1_mat'*tildeV_1_i;
            Z2(:,i)=Temp2_mat'*tildeV_1_i;
            Z3(:,i)=reshape(Temp1(:,i),[r,r])*tildeU_1_i;
            Z4(:,i)=reshape(Temp3(:,i),[n-r,r])*tildeU_1_i;
        end   
 
        B0=[Z1;Z2;Z3;Z4];

        % Bparas.B1=B1;
        % Bparas.B2=B2;
        % Bparas.B3=B3;
        % paras.m = m;
        % paras.n = n;
        % paras.r = r;

        %% cg method
        %runhist = cg_subproblem(Aparas,cg_paras,paras,B0,gamma);
        runhist = cg_subproblem_opt(Aparas,cg_paras,paras,B0,gamma);

        %fprintf("%d-th variable projection method, X1_Y1_cg iter=%d, X2_cg iter=%d, Y2_cg iter=%d\n",iter, runhist_X1_Y1.iter,runhist_X2.iter,runhist_Y2.iter);

        X = runhist.X; % [X',Y']'
        % X_opt=runhist_opt.X;
        % X_res=X_opt-X;
        % fprintf("X residual is %3.60f\n",norm(X_res));

        % A=PX 
        Atemp=X(1:m,1:r);
        A=P*Atemp;
        % gauss_tol=1e-5;
        % if norm(A(:))<gauss_tol
        %     A=10*A;
        % end 

        %fprintf("A norm is %3.10f\n",norm(A(:)));

        % B=QY   
        Btemp=X(m+1:m+n,1:r);
        B=Q*Btemp;
        %fprintf("B norm is %3.10f\n",norm(B(:)));
        % if norm(B(:))<gauss_tol
        %     B=10*B;
        % end 

        % % validate equation (17)
        % RRtemp=proj_supp(U,V)*(khatrirao(U,B)+khatrirao(A,V))*sym_inv(khatrirao(U,V));
        % LL=RRtemp*MA_temp-proj_supp(U,V)*MA_temp;
        % RR=RRtemp'*MA_temp;
        % result=trace(LL'*RR);
        % fprintf("============== result is %3.50f ============\n",result);

       %% check whether [A;B] is the descent direction
        dir=[A;B];
        % gra_gra_inner=-norm(gra(:))^2;
        % gra_direction_inner=trace(gra'*dir);
        % fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
    
     
        % if norm(dir(:))<epsilon
        %     fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
        %     break;
        % end    
    
        % if gra_direction_inner < gra_gra_inner 
        %     fprintf("-------------------------------- gauss-newton direction --------------------------------------------\n");
        % end    
            
        U_pre = U;
        V_pre = V;  
        
        fprintf("fval(U_pre,V_pre,MA_temp) is %3.8f\n",fval_pre);
        
        flag = false;
        alpha=1e-4;
        t=1.25;
        beta=0.8;
        threshold=alpha*t*trace(gra'*dir);
       
        for i=1:150
            t=beta*t;

            U=U_pre+t*A;
            V=V_pre+t*B;

            % update the fval
            UVkr=khatrirao(U,V);
            Pro=proj_supp(UVkr);
            PMA=Pro*MA_temp;
            
            %fval_cur=fval(U,V,MA_temp);
            fval_cur=fval_opt(PMA);


            [P, Sigma, S] = qr(U);
            S = S';
            P_tru=P(:,1:r);
            Sigma_tru=Sigma(1:r,:);

            [Q, Tau, T] = qr(V);
            T = T';

            Q_tru=Q(:,1:r);
            Tau_tru=Tau(1:r,:);
            X_hat=Sigma_tru*S;
            Y_hat=Tau_tru*T;
            XYkr=khatrirao(X_hat,Y_hat);

            [P_hat,Sigma_hat,Q_hat]=svd(XYkr,"econ");

            PQkr=kron(P_tru,Q_tru);
            PQP=PQkr*P_hat;
            Proj=PQP*PQP';

            % update the fval and the gradient
            R=MA_temp-Proj*MA_temp;
            fval_cur_opt=0.5*norm(R(:))^2;
    
            fval_cur_res=fval_cur-fval_cur_opt;
            fprintf("------------ the residual of fval_cur is %3.30f-------------\n",fval_cur_res);
    
   
 
            fprintf("i is %d, fval(U,V,MA)-fval(U_pre,V_pre,MA) is %3.8f, alpha*t*norm(gra(:))^2 is %3.8f,t is %3.8f\n", i, fval_cur-fval_pre,threshold,t);

            % if fval_cur-fval_pre<=0
            %     fprintf("================ fval_cur-fval_pre less than 0 and fval_cur-fval_pre is %3.20f\n",fval_cur-fval_pre);
            % end   
            
            % if fval_pre-fval_cur>0
            %     flag=true;
            %     fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
            %     break;
            % end    
            
            if fval_cur <= fval_pre+threshold
            %if fval_pre-fval_cur>0
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
                break;
            end
        end

        relfval=(fval_pre-fval_cur)/fval_init;
        UV=[U;V];
        relstep=t*norm(dir(:))/norm(UV(:));
        if relstep<1e-6
            fprintf("----------------------------- relstep is less than 1e-6 -------------------------\n");
            break;
        end
            
        fprintf("----------------------------- relfval is %3.20f, relstep is %3.20f -------------------------\n",relfval,relstep);
        if flag
            % update fval
            fval_pre=fval_cur;
            fval_pre_opt=fval_cur_opt;
            
            % update the gradient
            %gra=grad(U,V,MA_temp);
            gra=grad_opt(PMA,UVkr,U,V,MA_temp);

            UVi=Q_hat'*inv(Sigma_hat)*P_hat'*PQkr';
            PUVi=PQkr*P_hat*inv(Sigma_hat)*Q_hat;
            gra_opt_opt=grad_opt_opt_opt(R,UVi,PUVi,U,V,MA_temp);

            grad_res=gra_opt_opt-gra;
            fprintf("the residual of grad is %3.20f\n",norm(grad_res(:)));

        else
            U=U_pre;
            V=V_pre;
            fprintf("=================================line search failed!============================\n");
        end
        
        % fprintf("the norm of the gauss-newton direction is %3.10f, epsilon*(norm(Sol(:))+epsilon) is %3.10f\n",t*norm(dir(:)),epsilon*(norm(Sol(:))+epsilon));
        % if t*norm(dir(:))<epsilon*(norm(Sol(:))+epsilon)
        %     break;
        % end  

        fprintf("the norm of the gauss-newton direction is %3.10f\n",norm(dir(:)));
        if norm(dir(:))<epsilon
            fprintf("the algorithm is terminate by the gauss-newton direction condition\n");
            break;
        end

        if relfval <1e-12
            fprintf("--------------------------relfval less than 1e-12, the algorithm is terminate by the relfval condition---------------------\n");
            break;
        end  


        [P, Sigma, S] = qr(U);
        S = S';
        [Q, Tau, T] = qr(V);
        T = T';


        % if ~flag
        %     % gamma = 1;
        %     % U = U_pre + gamma * A;
        %     % V = V_pre + gamma * B;
        %     fprintf("line search failed!\n");
        % end

       
        
        fprintf("=========== the norm of the gradient is %3.8f ============\n",norm(gra(:)));

        fprintf("variable_projection_fopt %d-th iteration, fval is %3.20f =============================================================\n",iter,fval_cur);
        
    end

    krmp=pinv(UVkr); % the Moore-Penrose inverse of U\odot V
    W=MA_temp'*krmp';
    runhist.U = U;
    runhist.V = V;
    runhist.W = W;
    runhist.iter = iter;
    runhist.cput = cputime - timect;
end
