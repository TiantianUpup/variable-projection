function eqs = fval_tilde_opt_v3(WTW,N)
    NNT=N*N';
    
    eqs=-0.5*sum(WTW(:).*NNT(:));
end