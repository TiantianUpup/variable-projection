function H = hank(p,r)
    h = rand(p+r-1,1);
    H = hankel(h(1:p), h(p:end));
end    
