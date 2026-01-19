function runhist = mcg(A,b,paras)
    fprintf("-------------------------- Is A sparse %d ---------------\n",issparse(A));
    timect = cputime;
    itmax = paras.itmax;
    %tol=paras.tol;
    tol = paras.tol;

    % initialization
    [m,n]=size(b);

    fprintf("the size of b is\n");
    disp(size(b))
    x=zeros(m,n);
    fprintf("the size of x is\n");
    disp(size(x))

    r = b;

    iter = 0;
    r_norm = norm(r(:));

    while (iter < itmax && r_norm > tol)
        % step 1
        iter = iter + 1;

        % step 2
        if iter == 1
            c = r;
        else
            beta = norm(r(:))^2 / norm(r_pre(:))^2;
            c = r + beta * c;
        end

        % Step 3
        alpha = norm (r(:))^2 / trace(c'*A*c);

        % Step 4
        x = x + alpha*c;

        % Step 5
        r_pre = r;
        r = r-alpha*A*c;

        r_norm = norm(r(:));
    end

    cput = cputime - timect;
    %fprintf("cg method iter=%d, r_norm=%3.7f\n",iter,norm(r_norm(:)));
    runhist.r_norm = r_norm;
    fprintf("mcg method iter=%d, R_norm=%3.7f, time cost is %3.4f\n",iter,norm(r_norm),cput);
    runhist.x = x;
    runhist.iter = iter;


end