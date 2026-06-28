function runhist = variable_projection_simple(U0,V0,U_true,V_true,W_true,MA,vp_paras,cg_paras,paras)
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    [p,~]=size(W_true);
    [~,k] = size(U_true);
    fval=[];

    % U_0 and V_0
    U = U0;
    V = V0;
    WTW=W_true'*W_true;

    [P, Sigma, S] = qr(U);
    S = S';
    P1=P(:,1:r);
    [Q, Tau, T] = qr(V);
    T = T';
    Q1=Q(:,1:r);
    Sigma_1 = Sigma(1:r, 1:r);
    tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
    Tau_1 = Tau(1:r, 1:r);
    tildeV_1 = Tau_1 * T;
    G = khatrirao(tildeU_1, tildeV_1);
    [W, Db, F] = qr(G);
    W1=W(:,1:r);
    W2=W(:,r+1:end);
    W2W2T=W2*W2';
    F = F';
    D = Db(1:r, 1:r);
    PTU=P'*U_true;
    QTV=Q'*V_true;

    %gra=grad(U,V,MA);
    % UVkr=khatrirao(U,V);
    % Proj=proj_supp(UVkr);
    % PMA=Proj*MA;
    % gra=grad_opt(PMA,UVkr,U,V,MA);
    % grad_res=gra_opt-gra;
    % fprintf("the residual of grad is %3.20f\n",norm(grad_res(:)));

    %fval_pre=fval(U,V,MA);
    % fval_pre=fval_opt(PMA);
    % fval_init=fval_pre;

    % fval_res=fval_pre-fval_pre_opt;
    % fprintf("the residual of fval_pre is %3.30f\n",fval_res);

    % cput_finit=cputime;
    % UVkr=khatrirao(U,V);
    % Proj=proj_supp(UVkr);
    % PMA=Proj*MA;
    % fval_init=fval_opt(PMA);
    % fprintf("---------------------- the calculateion of the fval_init costs %3.10f ----------------\n",cputime-cput_finit);

    
    % cput_finit_opt=cputime;
    % tic
    % MA_norm=0;
    % for i=1:k
    %     Ui=U_true(:,i);
    %     Vi=V_true(:,i);
    %     Wi=W_true(:,i);
    %     for j=1:k
    %         Uj=U_true(:,j);
    %         Vj=V_true(:,j);
    %         Wj=W_true(:,j);
    %         MA_norm=MA_norm+(Ui'*Uj)*(Wi'*Wj)*(Vi'*Vj);
    %     end
    % end
    % toc
    % fprintf("---------------------- MA norm costs %3.10f ---------------\n",cputime-cput_finit_opt); 
    
    % cput_ma=cputime;
    % tic
    UTU=U_true'*U_true;
    VTV=V_true'*V_true;
    UTUVTV=UTU.*VTV;
    MA_norm=sum(UTUVTV(:).*WTW(:));
    % toc
    % fprintf("---------------------- MA_opt norm costs %3.10f ---------------\n",cputime-cput_ma); 

    % MA_res=MA_norm-MA_norm_opt;
    % fprintf("-------------------- MA_res is %3.20f -----------------\n",MA_res);
    
    fval_init=fval_tilde_opt(WTW,P1,Q1,U_true,V_true,W1)+0.5*MA_norm;
    %fprintf("---------------------- the calculateion of the fval_init costs %3.10f ----------------\n",cputime-cput_finit);
    fval_pre=fval_init;
    fprintf("variable_projection_simple method, the initial objective value fval is %3.4f\n",fval_init);

    iter = 0; 
    %itmax = paras.itmax;

    gamma=vp_paras.gamma;
    % fprintf("-----------------------------------------------gamma is %d\n",gamma);
    itmax=vp_paras.itmax;
    tol=vp_paras.tol; 
    epsilon=1e-6;
    cg_paras.mu=vp_paras.mu;
    
    % calculate the permutation matrix and its blocks
    cput_per=cputime;
    H = permutation(m,n,r);
    H1 = H(1:m*n,1:r^2);
    H21 = H(1:m*n,r^2+1:n*r);
    H22 = H(1:m*n,n*r+1:(n+m)*r-r^2);
    %fprintf("------------------------- the permutation method costs %3.10f --------------------------\n",cputime-cput_per);

    %&& norm(gra(:))>epsilon
    %MA_temp=MA;
    while (iter < itmax)
        iter = iter + 1;

        fval(end+1)=fval_pre;
        

        % gra=grad(U,V,MA);
        % check whether (U\odot V)-\in K(U\odot V)
        %UVkr=khatrirao(U,V);

        %fprintf("the rank of the UVkr is %d\n",rank(UVkr));
      
        % calculate \overline{M(A)}
        % Method 1
        % tic
        % kro = kron(P, Q);
        % MA = kro' * MA_temp;
        % toc
        % % Method 2
        % tic
        % MA_opt=zeros(n*m,p);
        % %QTV=Q'*V_true;
        % UTP=U_true'*P;
        % for i=1:p
        %     MA_opt(:,i)=reshape(QTV*diag(W_true(i,:))*UTP,[n*m,1]);
        % end    
        % toc
        % MA_res=MA-MA_opt;
        % fprintf("----------- the residual of MA is %3.40f ------------------------\n",norm(MA_res(:)));

        % calculate M
        HWJ = H1 * W1;
        Di=inv(D);% the dimension of D is r\times r
        M = F'*Di*HWJ';

        Beta=zeros(r,k);
        for i=1:r
            Beta(i,:)=sum(QTV.*(reshape(M(i,:),[n,m])*PTU),1);
        end

        PQkr = khatrirao(PTU,QTV)*WTW;
        MAC=PQkr*Beta';

        WWTH=W2W2T*H1';

        Temp1=WWTH*MAC;
        Temp2=H22'*MAC;
        Temp3=H21'*MAC;

        % calculate B
        Z1=zeros(r,r);
        Z2=zeros(m-r,r);
        Z3=zeros(r,r);
        Z4=zeros(n-r,r);

        for i=1:r
            tildeU_1_i=tildeU_1(:,i);
            tildeV_1_i=tildeV_1(:,i);

            Temp1_mat=reshape(Temp1(:,i),[r,r]);
            Temp2_mat=reshape(Temp2(:,i),[r,m-r]);

            Z1(:,i)=Temp1_mat'*tildeV_1_i;
            Z2(:,i)=Temp2_mat'*tildeV_1_i;
            Z3(:,i)=Temp1_mat*tildeU_1_i;
            Z4(:,i)=reshape(Temp3(:,i),[n-r,r])*tildeU_1_i;
        end   
 
        B0=[Z1;Z2;Z3;Z4];
        % B_res=B0-B0_simple;
        % fprintf("========================== the residual of B is %3.20f =======================\n",norm(B_res)); 

        Aparas.M=M;
        Aparas.MAC=MAC;
        Aparas.W2W2T=W2W2T;
        Aparas.tildeU_1=tildeU_1;
        Aparas.tildeV_1=tildeV_1;

        % add gamma term
        if gamma~=0
            % fprintf("----------------------------------------- gamma is %d\n",gamma);
            PQkr_1=kron(P,Q);
            MA_tilde=PQkr_1'*MA;

            Aparas.B1=W2'*H1'*MA_tilde;
            Aparas.B2=H21'*MA_tilde;
            Aparas.B3=H22'*MA_tilde;

            Aparas.E=Di'*F;
            Aparas.W2=W2;
        end

        runhist = cg_subproblem_simple(Aparas,cg_paras,paras,B0,gamma);

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
        dir=[A;B];
        %fprintf("B norm is %3.10f\n",norm(B(:)));
        % if norm(B(:))<gauss_tol
        %     B=10*B;
        % end 

        % projction A B onto the tangent space of Oblique manifold at U V
        A=proj_tan(A,U);
        B=proj_tan(B,V);

        % % validate equation (17)
        % RRtemp=proj_supp(U,V)*(khatrirao(U,B)+khatrirao(A,V))*sym_inv(khatrirao(U,V));
        % LL=RRtemp*MA_temp-proj_supp(U,V)*MA_temp;
        % RR=RRtemp'*MA_temp;
        % result=trace(LL'*RR);
        % fprintf("============== result is %3.50f ============\n",result);

        % check whether [A;B] is the descent direction
        
        % gra_gra_inner=-norm(gra(:))^2;
        % gra_direction_inner=trace(gra'*dir);
        % fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
    
     
        % if norm(dir(:))<epsilon
        %     fprintf("========================= gra_direction_inner is %3.10f, gra_gra_inner is %3.10f\n",gra_direction_inner,gra_gra_inner);
        %     break;
        % end    
    
        % fprintf("-------------------------------- gauss-newton direction norm is %3.10f, gradient direction, norm is %3.10f --------------------------------------------\n",norm(dir(:)),norm(gra(:)));
               
            
        U_pre = U;
        V_pre = V;  
        % P1_pre=P1;
        % Q1_pre=Q1;
        % W1_pre=W1;
                
        % flag = false;
        % %alpha=1e-4;
        % t=1.25;
        % beta=0.8;
        %threshold=alpha*t*trace(gra'*dir);


        % % test line search for the gradient method
        % gra_A=-gra(1:m,1:r);
        % gra_B=-gra(m+1:end,1:r);

        % for i=1:150
        %     t=beta*t;
        %     U=U_pre+t*gra_A;
        %     V=V_pre+t*gra_B;

        %     [P, Sigma, S] = qr(U);
        %     S = S';
        %     P1=P(:,1:r);
        %     [Q, Tau, T] = qr(V);
        %     T = T';
        %     Q1=Q(:,1:r);
        %     Sigma_1 = Sigma(1:r, 1:r);
        %     tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
        %     Tau_1 = Tau(1:r, 1:r);
        %     tildeV_1 = Tau_1 * T;
        %     G = khatrirao(tildeU_1, tildeV_1);
        %     [W, Db, F] = qr(G);
        %     W1=W(:,1:r);
        %     W2=W(:,r+1:end);
        %     W2W2T=W2*W2';
        %     F = F';
        %     D = Db(1:r, 1:r);

        %     PTU=P'*U_true;
        %     QTV=Q'*V_true;

        %     fval_cur=fval_tilde_opt(WTW,P1,Q1,U_true,V_true,W1);
        %     f_diff=fval_pre-fval_cur;

        %     if f_diff>0
        %         %if fval_pre-fval_cur>0
        %         %flag=true;
        %         fval_pre=fval_cur;
        %         fprintf("iter=%d, grad backtracking stepsize selection successful, step length is %3.8f,fval_cur is %3.10f\n",iter,t,fval_cur+0.5*MA_norm);
        %         break;
        %     % else
        %     %     fprintf("iter=%d-%d, grad line search failed\n",iter,i);
        %     end
        % end    

       
        %fval_pre = fval_tilde_opt(WTW,P1_pre,Q1_pre,U_true,V_true,W1_pre);
        t=1.25;
        beta=0.8;
        for i=1:150
            t=beta*t;

            U=U_pre+t*A;
            V=V_pre+t*B;

            [P, Sigma, S] = qr(U);
            S = S';
            P1=P(:,1:r);
            [Q, Tau, T] = qr(V);
            T = T';
            Q1=Q(:,1:r);
            Sigma_1 = Sigma(1:r, 1:r);
            tildeU_1 = Sigma_1 * S; % tildeU_1: \tilde{U}_1
            Tau_1 = Tau(1:r, 1:r);
            tildeV_1 = Tau_1 * T;
            G = khatrirao(tildeU_1, tildeV_1);
            [W, Db, F] = qr(G);
            W1=W(:,1:r);
            W2=W(:,r+1:end);
            W2W2T=W2*W2';
            F = F';
            D = Db(1:r, 1:r);

            PTU=P'*U_true;
            QTV=Q'*V_true;

            fval_cur=fval_tilde_opt(WTW,P1,Q1,U_true,V_true,W1);
            f_diff=fval_pre-fval_cur;
        
            %tic
            %fval_cur=fval_tilde_opt(WTW,P1,Q1,U_true,V_true,W1);
            % toc
            % tic
            % fdiff_opt=fval_tilde(P1,Q1,W1,MA_temp)-fval_tilde(P1_pre,Q1_pre,W1_pre,MA_temp);
            % toc
            % tic
            % fdiff=fval(U,V,MA_temp)-fval(U_pre,V_pre,MA_temp);
            % toc
            % f_gap=fdiff_opt-fdiff;
            % fprintf("---------------------------------- fdiff_opt_opt=%3.10f, fdiff_opt=%3.10f, fdiff=%3.10f --------------------\n",fdiff_opt_opt,fdiff_opt,fdiff);
            % fprintf("----------------------------------- f_gap=%3.40f ------------------------------------\n",f_gap);

            % UVkr=khatrirao(U,V);
            % Proj=proj_supp(UVkr);
            % PMA=Proj*MA_temp;
            
            %fval_cur=fval(U,V,MA_temp);
            % fval_cur=fval_opt(PMA);

            % fval_cur_res=fval_cur-fval_cur_opt;
            % fprintf("the residual of fval_cur is %3.30f\n",fval_cur_res);
    
   
 
            %fprintf("i is %d, fval(U,V,MA)-fval(U_pre,V_pre,MA) is %3.8f, alpha*t*norm(gra(:))^2 is %3.8f,t is %3.8f\n", i, fval_cur-fval_pre,threshold,t);

            % if fval_cur-fval_pre<=0
            %     fprintf("================ fval_cur-fval_pre less than 0 and fval_cur-fval_pre is %3.20f\n",fval_cur-fval_pre);
            % end   
            
            % if fval_pre-fval_cur>0
            %     flag=true;
            %     fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
            %     break;
            % end    
            
            if f_diff>0
                %if fval_pre-fval_cur>0
                %flag=true;
                fval_pre=fval_cur;
                %fprintf("iter=%d, vp backtracking stepsize selection successful, step length is %3.8f,fval_cur is %3.10f\n",iter,t,fval_cur+0.5*MA_norm);
                break;
            % else
            %    fprintf("iter=%d-%d, vp line search failed\n",iter,i);
            end
        end

        % normalize 
        U=proj_oblique(U);
        V=proj_oblique(V);

        relfval=f_diff/fval_init;
        UV=[U;V];
        relstep=t*norm(dir(:))/norm(UV(:));
        if relstep<=1e-6
            fprintf("----------------------------- relstep is less than 1e-6 -------------------------\n");
            break;
        end
            
        %fprintf("----------------------------- relfval is %3.20f, relstep is %3.20f -------------------------\n",relfval,relstep);
        % if flag
            
        %     % update fval
            
            
        %     %fval_pre=fval_cur;
            
        %     % update the gradient
        %     %gra=grad(U,V,MA_temp);
        %     %gra=grad_opt(PMA,UVkr,U,V,MA_temp);
        %     %grad_res=gra_opt-gra;
        %     %fprintf("the residual of grad is %3.20f\n",norm(grad_res(:)));

        % else
        %     U=U_pre;
        %     V=V_pre;
        %     fprintf("=================================line search failed!============================\n");
        % end
        
        % fprintf("the norm of the gauss-newton direction is %3.10f, epsilon*(norm(Sol(:))+epsilon) is %3.10f\n",t*norm(dir(:)),epsilon*(norm(Sol(:))+epsilon));
        % if t*norm(dir(:))<epsilon*(norm(Sol(:))+epsilon)
        %     break;
        % end  

        %fprintf("the norm of the gauss-newton direction is %3.10f\n",norm(dir(:)));
        % if norm(dir(:))<epsilon
        %     fprintf("the algorithm is terminate by the gauss-newton direction condition\n");
        %     break;
        % end

        if relfval <=1e-12
            fprintf("--------------------------relfval less than 1e-12, the algorithm is terminate by the relfval condition---------------------\n");
            break;
        end  



        % if ~flag
        %     % gamma = 1;
        %     % U = U_pre + gamma * A;
        %     % V = V_pre + gamma * B;
        %     fprintf("line search failed!\n");
        % end

       
        
        % fprintf("=========== the norm of the gradient is %3.8f ============\n",norm(gra(:)));

       % fprintf("variable_projection_gamma %d-th iteration, fval is %3.10f =============================================================\n",iter,fval_cur+0.5*MA_norm);
        
    end

    % gra=grad(U,V,MA);
    % fprintf("----------------------------- the gradient is %3.10f -----------------------------\n",norm(gra(:)));
    % fprintf("------------- the max norm of the gradient is %3.10f -----------------------------\n",max(abs(gra), [], 'all'));
    fval=fval+0.5*MA_norm;
    % tic
    % UVkr=khatrirao(U,V);
    % %cput_w=cputime;
    % krmp=pinv(UVkr); % the Moore-Penrose inverse of U\odot V
    % W=MA'*krmp';
    % toc
    % tic
    W=W_true*((U_true'*U).*(V_true'*V))*pinv((U'*U).*(V'*V));
    %toc
    %W_res=W-W_opt;
    %fprintf("------------------------ the residual of W is %3.20f -----------------------\n",norm(W_res(:)));
    %fprintf("--------------------------- calculate the W costs %3.10f -----------------------\n",cputime-cput_w);
    runhist.U = U;
    runhist.V = V;
    runhist.W = W;
    runhist.fval=fval;
    runhist.iter = iter;
    runhist.cput = cputime - timect;
end
