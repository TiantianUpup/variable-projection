function runhist = vp_pGN_ori(U0,V0,U_true,V_true,W_true,vp_paras,cg_paras)
                            
    timect = cputime;
    [m,r] = size(U0);
    [n,~] = size(V0);
    
    paras.m=m;
    paras.n=n;
    paras.r=r;
    
    % Running history
    fval=[];

    U = U0;
    V = V0;
    WTW=W_true'*W_true;
   
    [P, Sigma, S] = qr(U);
    S = S';
    P1=P(:,1:r);
    P2=P(:,r+1:m);
    [Q, Tau, T] = qr(V);
    T = T';
    Q1=Q(:,1:r);
    Q2=Q(:,r+1:n);
    Sigma_1 = Sigma(1:r, 1:r);
    tildeU_1 = Sigma_1 * S; 
    Tau_1 = Tau(1:r, 1:r);
    tildeV_1 = Tau_1 * T;
    G = khatrirao(tildeU_1, tildeV_1);
    [W, Db, F] = qr(G);
    W1=W(:,1:r);
    
    F = F';
    D = Db(1:r, 1:r);

    P1TU=P1'*U_true;
    P2TU=P2'*U_true;
    Q1TV=Q1'*V_true;
    Q2TV=Q2'*V_true;

    UTU=U_true'*U_true;
    VTV=V_true'*V_true;
    UTUVTV=UTU.*VTV;
    MA_norm=sum(UTUVTV(:).*WTW(:));
   
    N=cal_N(P1TU,Q1TV,W1);
    fval_init=fval_tilde(WTW,N)+0.5*MA_norm;
    fval_pre=fval_init;
    
    iter = 0;
    itmax=vp_paras.itmax;
    xtol=vp_paras.xtol;
    ftol=vp_paras.ftol; 
    
    while (iter < itmax)
        iter = iter + 1;
        fval(end+1)=fval_pre;

        Di=inv(D);
       
        % calculate B
        Z1=zeros(r,r);
        Z2=zeros(m-r,r);
        Z3=zeros(r,r);
        Z4=zeros(n-r,r);

        NDiF=N*(Di'*F);
        WTWNDiF=WTW*NDiF;
        CCT=NDiF'*WTWNDiF;

        for i=1:r
            tildeU_1_i=tildeU_1(:,i);
            tildeV_1_i=tildeV_1(:,i);
            WTWNDiFi=WTWNDiF(:,i);

            
            a=reshape(Q1TV.*WTWNDiFi'*P1TU',[r*r,1]);
            W2W2TH1TMACi=a-W1*(W1'*a);
            a_mat=reshape(W2W2TH1TMACi,[r,r]);

            Z1(:,i)=a_mat'*tildeV_1_i;
            Z2(:,i)=P2TU*(WTWNDiFi.*(Q1TV'*tildeV_1_i));
            Z3(:,i)=a_mat*tildeU_1_i;
            Z4(:,i)=Q2TV*(WTWNDiFi.*(P1TU'*tildeU_1_i));
        end   
 
        B1=[Z1;Z3];

        Aparas.NDiF=NDiF;
        Aparas.WTWNDiF=WTWNDiF;
        Aparas.W1=W1;
        Aparas.tildeU_1=tildeU_1;
        Aparas.tildeV_1=tildeV_1;
        Aparas.CCT=CCT;

        runhist_X1_Y1 = X1_Y1_subproblem(Aparas, B1, cg_paras, paras);
        runhist_X2 = X2_subproblem(Aparas, Z2, cg_paras, paras);
        runhist_Y2 = Y2_subproblem(Aparas, Z4, cg_paras, paras);

        X1_Y1 = runhist_X1_Y1.X; % [X',Y']'
        X2 = runhist_X2.X2;
        Y2 = runhist_Y2.Y2;
        X1=X1_Y1(1:r,1:r);
        Y1=X1_Y1(r+1:2*r,1:r);
        Atemp=[X1;X2];
        Btemp=[Y1;Y2];

        A=P*Atemp;
        B=Q*Btemp;
        dir=[A;B];
        
        U_pre = U;
        V_pre = V;  
        
       
        t=1.25;
        beta=0.8;
        for i=1:150
            t=beta*t;

            U=U_pre+t*A;
            V=V_pre+t*B;

            [P, Sigma, S] = qr(U);
            S = S';
            P1=P(:,1:r);
            P2=P(:,r+1:m);
            [Q, Tau, T] = qr(V);
            T = T';
            Q1=Q(:,1:r);
            Q2=Q(:,r+1:n);
            Sigma_1 = Sigma(1:r, 1:r);
            tildeU_1 = Sigma_1 * S; 
            Tau_1 = Tau(1:r, 1:r);
            tildeV_1 = Tau_1 * T;
            G = khatrirao(tildeU_1, tildeV_1);
            [W, Db, F] = qr(G);
            W1=W(:,1:r);
            F = F';
            D = Db(1:r, 1:r);

            P1TU=P1'*U_true;
            P2TU=P2'*U_true;
            Q1TV=Q1'*V_true;
            Q2TV=Q2'*V_true;

            N=cal_N(P1TU,Q1TV,W1);

            fval_cur=fval_tilde(WTW,N);
            f_diff=fval_pre-fval_cur;
        
           
            if f_diff>0
                fval_pre=fval_cur;
                break;
            end
        end

        % Projction U, V onto Oblique manifold.
        U=proj_oblique(U);
        V=proj_oblique(V);

        % Check for convergence.
        UV=[U;V];
        relstep=t*norm(dir(:))/norm(UV(:));
        if relstep<=xtol
            runhist.info=2;
            break;
        end
            
        relfval=f_diff/fval_init;
        if relfval <=ftol
            runhist.info=1;
            break;
        end  
    end

    if iter >= itmax
        runhist.info=1;
    end     

    fval=fval+0.5*MA_norm;
    W=W_true*((U_true'*U).*(V_true'*V))*pinv((U'*U).*(V'*V));
    

    % Update the output structure.
    runhist.U = U;
    runhist.V = V;
    runhist.W = W;
    runhist.fval=fval;
    runhist.iter = iter;
    runhist.cput = cputime - timect;
end
