addpath("E:\\matlab-code\\variable-projection\\tensor_toolbox")

U = [1,2;
     3,4];
V = [1,1;
     1,1;
     1,1];
W = [5,6;
     7,8];

% calculation of MA
kr = khatrirao(U, V);
disp(kr);
% MA = kr * W';
% disp(MA);

