function X_noise = add_laplace_noise(X,epsilon)

    N = laplace_noise(size(X,1),size(X,2),1);

    N = epsilon*norm(X,'fro')/norm(N,'fro')*N;

    X_noise = X + N;

end

function N = laplace_noise(m,n,b)
    u = rand(m,n)-0.5;
    N = -b * sign(u) .* log(1-2*abs(u));
end