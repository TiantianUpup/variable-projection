function P = perm_cg(m,n,r)
% 构造置换矩阵 P，使得 W = P * vec(Z)
% 输入:
%   r: X_1, Y_1 的维度 (r×r)
%   m: Z 的总行数的一部分 (X_1 和 X_2 的总行数)
%   n: Z 的总行数的另一部分 (Y_1 和 Y_2 的总行数)
% 输出:
%   P: 置换矩阵 (N×N), N = r*(m+n)

M = m + n;          % Z 的总行数
N = r * M;          % 总维度

% 初始化置换矩阵
P = zeros(N, N);

for q = 1:r         % 列索引 j = 1,...,r
    % 情况1: X_1 (行 1..r)
    for p = 1:r
        k = (q-1)*M + p;        % vec(Z) 中的索引 (1-based)
        w = (q-1)*r + p;        % W 中的索引 (1-based)
        P(w, k) = 1;
    end
    
    % 情况2: X_2 (行 r+1..m)
    for p = r+1:m
        k = (q-1)*M + p;
        w = r^2 + (q-1)*(m-r) + (p-r);
        P(w, k) = 1;
    end
    
    % 情况3: Y_1 (行 m+1..m+r)
    for p = m+1:m+r
        k = (q-1)*M + p;
        w = r^2 + r*(m-r) + (q-1)*r + (p-m);
        P(w, k) = 1;
    end
    
    % 情况4: Y_2 (行 m+r+1..m+n)
    for p = m+r+1:m+n
        k = (q-1)*M + p;
        w = 2*r^2 + r*(m-r) + (q-1)*(n-r) + (p-m-r);
        P(w, k) = 1;
    end
end

% % 验证 P 是置换矩阵
% assert(all(sum(P, 1) == 1), '每列应恰有一个1');
% assert(all(sum(P, 2) == 1), '每行应恰有一个1');
% assert(norm(P'*P - eye(N)) < 1e-10, 'P 应是正交矩阵');
end