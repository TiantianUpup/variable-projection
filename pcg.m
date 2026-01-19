function runhist = pcg(A,b,paras)
    timect = cputime;
    itmax = paras.itmax;
    %tol=paras.tol;
    tol = paras.tol;

    % initialization
    [m,n]=size(b);
    x=zeros(m,n);

    r = -b;

    % Solve My=r, y=M^{-1}r
    A=sparse(A);
    L=ichol(A); % A=L*L'
    z=L\r;
    y=L'\z;

    p=-y;

    iter = 0;
    r_norm = norm(r(:));

    while (iter < itmax && r_norm > tol)
       
        iter = iter + 1;

        alpha = (r'*y) / trace(p'*A*p);
        x=x+alpha*p;
        r_pre=r;
        r=r+alpha*A*p;

        % Solve My=r, y=M^{-1}r
        y_pre=y;
        z=L\r;
        y=L'\z;

        beta=-(r'*y)/(r_pre'*y_pre); 
        p=-y+beta*p;

        r_norm = norm(r(:));
    end

    cput = cputime - timect;
    %fprintf("cg method iter=%d, r_norm=%3.7f\n",iter,norm(r_norm(:)));
    runhist.r_norm = r_norm;
    fprintf("pcg method iter=%d, R_norm=%3.7f, time cost is %3.4f\n",iter,norm(r_norm),cput);
    runhist.x = x;
    runhist.iter = iter;


end