% data = load('T.mat');
% T = data.T;
% T_norm = normalize(T, 2,'range',[0,1]); 
% save('T-norm.mat', 'T_norm');


% data = load('T-norm.mat');
% T = data.T_norm;
% disp(T);


% data = load('R.mat');
% R = data.R;
% disp(R);

data = load('R.mat');
R = data.R;
R_norm = normalize(R, 2,'range',[0,1]); 
disp(R_norm);
save('R-norm.mat', 'R_norm');