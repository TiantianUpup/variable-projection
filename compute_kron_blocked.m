function M = compute_kron_blocked(P, V)
    % [P_rows, m] = size(P);
    % [r, r_cols] = size(V);
    
    % % 快速验证
    % if r_cols ~= r || P_rows ~= r*r
    %     error('维度错误: P应为%d×%d, V应为%d×%d', r*r, m, r, r);
    % end
    
    % % 预分配结果
    % M = zeros(r*r*r, m);
    
    % % 计算每个块
    % for i = 1:r
    %     % 提取 P_i
    %     row_start = (i-1)*r + 1;
    %     row_end = i*r;
    %     P_i = P(row_start:row_end, :);
        
    %     % 提取 v_i
    %     v_i = V(:, i);
        
    %     % 计算块并放入结果
    %     block_start = (i-1)*(r*r) + 1;
    %     block_end = i*(r*r);
    %     M(block_start:block_end, :) = kron(P_i, v_i);
    % end
end