
W=[1,2;
   3,4;
   5,6;
   7,8];
r=2;   
W_block = mat2cell(W, repmat(r, r, 1), r);

disp(W_block{1});
disp("============================");
disp(W_block{2});