%% this function calculates C^*C(B)
function X = BadjointB (tildeU_1, MMA, paras, C)
    
    n=paras.n;
    r=paras.r;
    X = zeros(n-r,r);
    I = eye(n-r);
    MMAt = MMA'; % the transpose of the MMA
    kr=khatrirao(tildeU_1,C);
    for i=1:r
        tildeU_1_i=tildeU_1(:,i);
        X(:,i)=kron(tildeU_1_i',I)*kr*MMA*MMAt(:,i);
    end 
end     