%function [x, i, res] = conjugate_gradient(A, b, eta, i_max)
% CONJUGATE_GRADIENT 共轭梯度法求解 Ax = b
% 输入:
%   A      - 对称正定矩阵
%   b      - 右端向量
%   eta    - 残差容限，当 ||r|| <= eta 时停止
%   i_max  - 最大迭代次数
% 输出:
%   x      - 近似解
%   i      - 实际迭代次数
%   res    - 每次迭代的残差范数记录

% 初始化
x = zeros(size(b));
r = b;
i = 0;
res = [norm(r)];

% 主循环
while (norm(r) > eta) && (i < i_max)
    i = i + 1;
    
    % Step 1.2: 计算搜索方向 p
    if i == 1
        p = r;
    else
        beta = (r_prev' * r_prev) / (r_prev_prev' * r_prev_prev);
        p = r + beta * p_prev;
    end
    
    % Step 1.3: 计算步长 alpha
    Ap = A * p;
    alpha = (r' * r) / (p' * Ap);
    
    % Step 1.4: 更新解
    x = x + alpha * p;
    
    % Step 1.5: 更新残差
    r_prev_prev = r_prev;  % 保存前两次的残差
    r_prev = r;
    r = r - alpha * Ap;
    
    % 保存当前迭代的残差范数
    res = [res; norm(r)];
    
    % 保存上一次的搜索方向
    p_prev = p;
end

%end