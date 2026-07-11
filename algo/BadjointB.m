%% this function calculates C^*C(B)
function X = BadjointB (Aparas,paras,J,gamma)
    CCT=Aparas.CCT;
    tildeV_1=Aparas.tildeV_1;

    m=paras.m;
    n=paras.n;
    r=paras.r;
   
    
    J2=J;
   
    for i=1:r
        CCi=CCT(:,i);
        
        tildeV_1_i=tildeV_1(:,i);
        X(:,i)=J2*(CCi.*(tildeV_1'*tildeV_1_i));
        
    end    

   
end     