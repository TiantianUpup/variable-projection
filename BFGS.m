function runhist = BFGS(U0,V0,MA, paras)
    fprintf("this is BFGS method\n");
    timect = cputime;
    m = paras.m;
    n = paras.n;
    r = paras.r;
    % U_0 and V_0
    U = U0;
    V = V0;
    gra=grad(U,V,MA);
    A=-gra(1:m,1:r);
    B=-gra(m+1:m+n,1:r);

    U_pre=U;
    V_pre=V;
    alpha=1e-4;
    t=1.25;
    beta=0.8;
    for i=1:150
        t=beta*t;

        U=U_pre+t*A;
        V=V_pre+t*B;
        
        if fval(U,V,MA)<=fval(U_pre,V_pre,MA)-alpha*t*norm(gra(:))^2
            flag=true;
            fprintf("backtracking stepsize selection successful, step length is %3.8f\n",t);
            break;
        end
    end

    gra_cur=grad(U,V,MA);
    S=t*[A;B];
    Y=gra_cur-gra;

    I=eye(m+n);
    H=trace(Y'*S)/norm(Y(:))^2*I;

    fprintf("beta is %3.4f\n",trace(Y'*S)/norm(Y(:))^2);
    fprintf("this is BFGS method, the initial objective value fval is %3.4f\n",fval(U,V,MA));

    % calculate B
    %B=-grad(U,V,MA);
    
    iter=0;

    itmax=1500;
    epsilon=1e-6;
    while (iter < itmax && norm(gra(:))>epsilon)
        iter=iter+1;
        % compute the search direction
        P=-H*gra;
        A=P(1:m,1:r);
        B=P(m+1:m+n,1:r);

        gra_direction_inner=trace(gra'*P);
        fprintf("the inner product of the direction and gradient is %3.4f\n",gra_direction_inner);

        %gra_pre=gra;

        % line search procedure (Wolfe condition)
        U_pre=U;
        V_pre=V;
        alpha=1.25;
        beta=0.8;
        c1=1e-4;
        c2=0.9;
        for i=1:150
            alpha=beta*alpha;

            U=U_pre+alpha*A;
            V=V_pre+alpha*B;

            gra_cur=grad(U,V,MA);
           
            fprintf("i is %d, fval(U,V,MA)-fval(U_pre,V_pre,MA) is %3.8f, c1*alpha*trace(gra'*P) is %3.8f,t is %3.8f\n", i, fval(U,V,MA)-fval(U_pre,V_pre,MA), c1*alpha*trace(gra'*P),alpha);
            fprintf("i is %d, abs(trace(gra_cur'*P)) is %3.8f, c2*abs(trace(gra'*P) is %3.8f\n", i,trace(gra_cur'*P), c2*trace(gra'*P));
            % fprintf("i is %d, fval(U_pre,V_pre,MA)-fval(U,V,MA) is %3.8f, alpha*t*norm(gra(:))^2 is %3.8f,t is %3.8f\n", i, fval(U_pre,V_pre,MA)-fval(U,V,MA), alpha*t*norm(gra(:))^2,t);
            %fprintf("trace(gra_cur'*P) is %3.4f, c2*trace(gra'*P) is %3.4f\n",trace(gra_cur'*P),c2*trace(gra'*P));   
            if fval(U,V,MA)-fval(U_pre,V_pre,MA)<=c1*alpha*trace(gra'*P)
                %&& abs(trace(gra_cur'*P))<=c2*abs(trace(gra'*P))
                % && trace(gra'*P)>=c2*gra_direction_inner
                flag=true;
                fprintf("backtracking stepsize selection successful, step length is %3.8f\n",alpha);
                break;
            end
        end

        gra_pre=gra;
        gra=grad(U,V,MA);

        % update the H
        S=t*[A;B];
        Y=gra-gra_pre;
        rho=1/trace(Y'*S);
        H=(I-rho*S*Y')*H*(I-rho*Y*S')+rho*(S*S');
    end

    runhist.iter=iter;
    runhist.U=U;
    runhist.V=V;
    runhist.cput=cputime-timect;
end    