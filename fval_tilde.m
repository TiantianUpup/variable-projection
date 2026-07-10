function eqs = fval(WTW,N)
    NNT=N*N';
    
    eqs=-0.5*sum(WTW(:).*NNT(:));
end