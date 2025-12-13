%% this function calculates C^*C(B)
function X = CadjointC (tildeV_1, MMA, paras, C)
    m=paras.m;
    r=paras.r;
    X = zeros(m-r,r);
    I = eye(m-r);
    MMAt = MMA'; % the transpose of the MMA
    kr=khatrirao(C,tildeV_1);
    for i=1:r
        tildeV_1_i=tildeV_1(:,i);
        X(:,i)=kron(I,tildeV_1_i')*kr*MMA*MMAt(:,i);
    end 
end    