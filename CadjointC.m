%% this function calculates C^*C(B)
function X = CadjointC (Aparas,paras,J,gamma)
    CCT=Aparas.CCT;

    tildeU_1=Aparas.tildeU_1;

    m=paras.m;
    n=paras.n;
    r=paras.r;
   
    J4=J;

    for i=1:r
        CCi=CCT(:,i);
        tildeU_1_i=tildeU_1(:,i);
     
        X(:,i)=J4*(CCi.*(tildeU_1'*tildeU_1_i)); 
    end    
end    