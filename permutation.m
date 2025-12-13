% this function generate a permutation matrix 
function H = permutation(m,n,r)
    I = eye(m*n,m*n);
   
    B = [];
    for k = 1 : r-1
        B = [B, k*n+1 : k*n+r];
    end

    C = [];
    for k = r : m-1
        C = [C, k*n+1 : k*n+r];
    end

    A = 1 : r;
    D1 = setdiff(1:n*r, [A, B]);
    D2 = setdiff(n*r+1:m*n, C);

    new_order = [A, B, D1, C, D2];
    H= I(:,new_order);
    %disp(new_order);

    % % moving columns
    % move_cols = [];
    % for k = 1 : r-1
    %     move_cols = [move_cols, k*n+1 : k*n+r];
    % end

    % % Not being moved and not in column 1:r
    % other_cols = setdiff(all_cols, [1:r, move_cols]);

    % % new indexs
    % new_order = [1:r, move_cols, other_cols];

    % disp(new_order);
    % H= I(:,new_order);

    % k = 1;
    % for i=1:(r-1)
    %     for j=1:r
    %         temp=H(:,r+k);
    %         H(:,r+k)=H(:,i*n+j);
    %         H(:,i*n+j)=temp;
 

    %         % temp=H(i*n+j,:);
    %         % H(i*n+j,:)=H(r+k,:);
    %         %H(r+k,:)=temp;
    %         k=k+1;
    %     end    
    % end    

    % k = 1;
    % for i=1:(r-1)
    %     for j=1:r
    %         temp=H(:,(r+i)*n+j);
    %         H(:,(r+i)*n+j)=H(:,r*(n+1)+k);
    %         H(:,r*(n+1)+k)=temp;

    %         k=k+1;
    %     end    
    % end    
end     