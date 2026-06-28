clc; close all;
% %%%%%%%% 100*100*100
% cp_als = [0.17, 2.35, 3.12, 5.12, 5.51,...
% 4.93, 5.67, 5.52, 5.32, 7.09, 7.26, 5.02, 6.50, 6.89, 13.55];
% cpd_als = [0.22, 1.65, 2.89, 3.64, 3.34,...
% 3.18, 3.50, 4.08, 3.44, 4.15, 4.88, 3.19, 3.46, 4.06, 4.42];
% cpd_nls = [0.42, 8.73, 8.63, 11.25, 15.51,...
% 10.60, 19.23, 12.88, 17.87, 12.34, 19.46, 12.83, 13.55, 12.45, 2.36];
% vp_pgn = [0.03, 0.73, 1.45, 2.03, 2.51,...
%  1.63, 3.21, 2.32, 3.18, 3.01, 5.75, 3.12, 2.66, 2.41, 0.90];

%%%%%%%% 30*40*1000
% cp_als = [0.22, 4.26, 6.11, 7.57, 10.66,...
% 11.44, 12.61, 13.05, 16.65, 15.90, 16.59, 18.10, 15.47, 19.44, 32.07];
% cpd_als = [0.18, 2.57, 3.13, 4.67, 5.56,...
% 5.58, 6.59, 6.74, 7.56, 7.12, 7.94, 8.53, 8.53, 8.85, 9.16];
% cpd_nls = [0.43, 7.41, 6.60, 9.32, 14.53,...
% 25.16, 28.04, 67.14, 103.35, 178.35, 134.31, 168.22, 102.82, 75.99, 26.50];
% vp_pgn = [0.03, 0.10, 0.14, 0.26, 0.43,...
% 0.73, 0.68, 1.10, 0.86, 1.17, 1.34, 1.54, 1.18, 0.85, 0.75];


%%%%%%%%%  100*100*10000
cp_als = [6.16, 160.60, 337.20, 543.08, 515.57,...
716.06, 878.89, 835.48, 863.08, 978.23, 924.95, 837.54, 1033.44, 933.45, 1338.52];
cpd_als = [5.54, 112.68, 217.15, 249.14, 261.16,...
280.26, 380.46, 356.30, 444.07, 532.40, 436.10, 408.93, 484.32, 441.88, 485.26]; 
cpd_nls = [17.22, 568.68, 575.53, 698.89, 1020.88,...
1402.46, 1726.92, 1414.34, 1202.03, 1540.14, 1242.34, 977.36, 1391.68, 982.34, 205.87];
vp_pgn= [0.05, 0.76, 1.35, 2.50, 3.86,...
4.55, 6.36, 4.00, 7.95, 8.41, 4.63, 3.50, 7.45, 4.81, 0.86]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot
fig=figure();
x=1:15;
% semilogy(x,cp_als,'-o','LineWidth', 1.5); hold all; %'r-o',
% semilogy(x,cpd_als,'-*', 'LineWidth', 1.5); hold all; %'g--s',
% semilogy(x,cpd_nls,'-square', 'LineWidth', 1.5); hold all; % 'g--s',
% semilogy(x,vp_pgn,'-x', 'LineWidth', 1.5);  % 'b-*',

plot(cp_als,'-o','LineWidth', 1.5); hold all; %'r-o',
plot(cpd_als,'-*', 'LineWidth', 1.5); hold all; %'g--s',
plot(cpd_nls,'-square', 'LineWidth', 1.5); hold all; % 'g--s',
plot(vp_pgn,'-x', 'LineWidth', 1.5);  % 'b-*',

% semilogy(fval_vp); hold all;
% semilogy(fval_cpd_als); hold all;
% semilogy(fval_cpd_nls);
ylabel('CPU time','FontSize',20,'FontWeight','bold'); 
xlabel('rank r','FontSize',20,'FontWeight','bold');
set(gca,'FontSize',13,'FontWeight','bold');
%title('Convergence plot'); 
legend('cp\_als','cpd\_als','cpd\_nls','vp\_pGN','FontSize', 15,'FontWeight', 'bold','Location', 'northwest');  % 

xlim([1, 15]);
xticks([1, 5, 10, 15]);
hold off;
% grid on
% grid minor off


ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - 1.1*ti(2) -1.1* ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,'p-10000','-dpdf') 
