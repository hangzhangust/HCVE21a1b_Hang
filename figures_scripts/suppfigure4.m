%% autocorrelation D=30
step =50;
run startup.m
load('Autocorrelation_new.mat', 'a_30')
load('Autocorrelation_new.mat', 'b_30')
% a_30 = cal_autocorrelation(all_Energy_1a_30,step);
% 
% b_30= cal_autocorrelation(all_Energy_1b_30,step);

FIG=figure;
p1=plot(a_30,'Color',purple,'LineWidth',1);
hold on
p2=plot(b_30,'Color',orange,'LineWidth',1)

set(gca,'XTick', [0 25 50])
xlim([0 50])
xlabel('Number of random walk steps, {\it k}')
ylabel({'Autocorrelation'})
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'autocorrelation';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 6.98 4]);
set(gca,'Position',[.18 .25 .78 .71]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
legendflex([p1,p2], {'Subtype 1a','Subtype 1b'}, 'ref', gcf, ...
                       'anchor', {'ne','ne'}, ...
                       'buffer',[-10 0], ...
                       'nrow',2, ...
                       'fontsize',8,'box','off','xscale',0.5);
set(gca,'TickLength',[0.02, 0.03])
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600'); end

%%
load('Neutrality_new.mat')
% epsilon=0:9;
% x = 10.^epsilon*10^-6;
% epsilon = 10^-5:10^-5:10^-3;
% epsilon =  [1e-02 2.5e-2 5e-02 7.5e-2 1e-01 2.5e-1 5e-1 7.5e-1 1 2.5 5 7.5 10];
% epsilon =  [1e-01 3e-01 5e-01 8e-01 3 10 30 50 80];
FIG=figure;
% loglog(epsilon,all_d_1a_L_500,'Color',purple);
% hold on
% loglog(epsilon,all_d_1b_L_500,'Color',orange)
p1=loglog(epsilon,all_d_1a_L_1000,'Color',purple,'LineWidth',1);
hold on
p2=loglog(epsilon,all_d_1b_L_1000,'Color',orange,'LineWidth',1)

% loglog(epsilon,all_d_1a_L_1000,'Color',purple);
% hold on
% loglog(epsilon,all_d_1b_L_1000,'Color',orange);


% loglog(epsilon,all_d_1a_L_1000,'LineStyle' ,'--','Color',purple);
% hold on
% loglog(epsilon,all_d_1b_L_1000,'LineStyle' ,'--','Color',orange);
% legend('Subtype 1a', 'Subtype 1b','Location','best')
xlabel('Acceptable change in fitness, \epsilon')
ylabel({'Neutrality'})
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'neutrality';
box off
FIG.Units = 'centimeters';
set(gca,'XTick', [0.01 1 10])
set(gca,'YTick', [0.01 1 100])
set(gcf,'Position',[10 10 6.98 4]);
set(gca,'Position',[.18 .25 .78 .71]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickDir','out')
% legendflex([p1,p2], {'Subtype 1a','Subtype 1b'}, 'ref', gcf, ...
%                        'anchor', {'ne','ne'}, ...
%                        'buffer',[-10 0], ...
%                        'nrow',2, ...
%                        'fontsize',8,'box','off','xscale',0.5);
set(gca,'TickLength',[0.02, 0.03])
set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.2, 0]);
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.18, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');  

% end
