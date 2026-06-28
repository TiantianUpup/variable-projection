function X = proj_tan(A,U)
    % the dimension of A and U is the same
    [m,r]=size(A);

    X=zeros(m,r);

    for i=1:r 
        Ai=A(:,i);
        Ui=U(:,i);

        X(:,i)=Ai-(Ui'*Ai)*Ui;
    end    
end