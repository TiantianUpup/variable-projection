%% this function calculates C^*C(B)
function X = CadjointC (Aparas,paras,J)
    CCT=Aparas.CCT;
    tildeU_1=Aparas.tildeU_1;
    r=paras.r;
   
    for i=1:r
        CCi=CCT(:,i);
        tildeU_1_i=tildeU_1(:,i);
     
        X(:,i)=J*(CCi.*(tildeU_1'*tildeU_1_i)); 
    end    
end    