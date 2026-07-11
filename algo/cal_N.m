% this function calculate N
function N=cal_N(P1TU,Q1TV,W1)
    [~,r]=size(W1);
    [~,q]=size(P1TU);

    N=zeros(q,r);
    
    for i=1:r
        wi=W1(:,i);
        W{i}=reshape(wi,[r,r]);

        for j=1:q
            P1TUj=P1TU(:,j);
            Q1TVj=Q1TV(:,j);
            N(j,i)=Q1TVj'*W{i}*P1TUj;
        end    
    end    
end    