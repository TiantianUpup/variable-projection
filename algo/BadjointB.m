%% this function calculates C^*C(B)
function X = BadjointB (Aparas,paras,J)
    CCT=Aparas.CCT;
    tildeV_1=Aparas.tildeV_1;
    r=paras.r;
   
    for i=1:r
        CCi=CCT(:,i);
        
        tildeV_1_i=tildeV_1(:,i);
        X(:,i)=J*(CCi.*(tildeV_1'*tildeV_1_i));     
    end    
end     