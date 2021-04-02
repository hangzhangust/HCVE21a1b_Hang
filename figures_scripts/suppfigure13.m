%%
load('mean_escape_time_1b.mat')
all_E=all_mean_escape_time_1b;

escape = [384;386;388;390;391;393;394;...
    395;396;397;398;399;400;401;402;403;404;405;407;408;410;415;416;417;422;424;431;433;...
    434;435;438;442;444;446;453;456;461;466;475;482;501;524;528;531;533;538;557;558;560;580;608;610;636;713];
escape = escape-383;
Escape = all_E(escape);
Remain= all_E;
Remain(escape)=[];

TPR=[];
FPR=[];
MCC=[];
F1=[];
for thresh = min(all_E):1:max(all_E)+1
TP = sum(Escape<thresh);
FN = sum(Escape>=thresh);
FP = sum(Remain<thresh);
TN = sum(Remain>=thresh);
precision = TP/(TP+FP); 
recall = TP/(TP+FN);
TPR = [TPR;TP/(TP+FN)];
FPR = [FPR;FP/(FP+TN)];
F1 = [F1;2*recall*precision/(recall+precision)];
MCC = [MCC;(TP*TN-FP*FN)/((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))^0.5];


end
FIG=figure
plot(FPR,TPR,'Color','k','LineWidth',1)

xlabel('False Positive Rate')
ylabel({'True Positive Rate'})
FIG.Units = 'centimeters';
FIG.Name = 'ROC';
set(gcf,'Position',[10 10 5 4]);

set(gca,'Position',[.2 .2 .74 .71]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
% text(0.6,0.1,sprintf('$$ AUC = %.2f$$',trapz(FPR,TPR)),'interpreter','latex','FontSize',8)

set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% figure;
% group = [zeros(size(Escape)),ones(size(Remain))];
% boxplot([Escape,Remain],group,'Labels',{'Escape','Remain'})
FIG=figure
[~,ind] = max(MCC);
[~,ind1] = max(F1);


x = min(all_E):1:max(all_E)+1;
plot(x,MCC,'Color','k','LineWidth',1)
hold on;
color = [0.6000    0.6000    0.6000];%gray
plot(x,F1,'Color',color,'LineWidth',1);
hold on;
yt = get(gca, 'YTick');
plot([x(ind) x(ind)], [min(yt) max(yt)], '-k','Color','k','LineWidth',1)

xlabel('Escape time cut-off')
ylabel({'True Positive Rate'})
FIG.Units = 'centimeters';
FIG.Name = 'F1_1b';
set(gcf,'Position',[10 10 5 4]);
set(gca,'Position',[.2 .2 .74 .71]);  %调整 XLABLE和YLABLE不会被切掉
ylim([0 0.8])
% text(x(ind)-20,max(yt)*1.05,sprintf('$$ %.2f$$',x(ind)),'interpreter','latex','FontSize',12)
legendflex(gca, {'MCC','F1 score'}, 'xscale', 0.5,'box','off');
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.2, 0.5, 0]);
% set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.5, -0.08, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end