%% Autocorrelation
step =50;
run startup.m
load('Autocorrelation_new.mat', 'a_5')
load('Autocorrelation_new.mat', 'b_5')
% a_5 = cal_autocorrelation(all_Energy_1a_5,step);
% b_5= cal_autocorrelation(all_Energy_1b_5,step);
FIG=figure;
p1=plot(a_5,'Color',purple,'LineWidth',1);
hold on
p2=plot(b_5,'Color',orange,'LineWidth',1)
% plot(a_30,'LineStyle' ,'--','Color',purple);
% hold on
% plot(b_30,'LineStyle' ,'--','Color',orange);
% legend('Subtype 1a', 'Subtype 1b','Location','best')
% legend boxoff

set(gca,'XTick', [0 25 50])
xlim([0 50])
xlabel('Number of random walk steps, {\it k}')
ylabel({'Autocorrelation'})
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'autocorrelation';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 6.98 5]);
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
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      end

%% 
load('Neutrality_new.mat')
% epsilon=0:9;
% x = 10.^epsilon*10^-6;
% epsilon = 10^-5:10^-5:10^-3;
% epsilon =  [1e-02 2.5e-2 5e-02 7.5e-2 1e-01 2.5e-1 5e-1 7.5e-1 1 2.5 5 7.5 10];
% epsilon =  [1e-01 3e-01 5e-01 8e-01 3 10 30 50 80];
FIG=figure;
p1=loglog(epsilon,all_d_1a_L_500,'Color',purple,'LineWidth',1);
hold on
p2=loglog(epsilon,all_d_1b_L_500,'Color',orange,'LineWidth',1)

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
set(gcf,'Position',[10 10 6.98 5]);
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
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 
% catch      
% 
% end

%% Peak analysis

%escape mutations
all_sites=[384;386;388;390;391;393;394;...
    395;396;397;398;399;400;401;402;403;404;405;407;408;410;415;416;417;422;424;431;433;...
    434;435;438;442;444;446;453;456;461;466;475;482;501;524;528;531;533;538;557;558;560;580;608;610;636;713];


load('Peaks.mat')
load('Model_1a.mat', 'num_mutants_combine_array');
load('Model_1a.mat', 'ind_non_conserve')
load('Model_1a.mat', 'msa_aa')
load('Model_1a.mat', 'J_mat')
%1a
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
J_mat = reshape(J_mat,[624 624]);
J_MPF_BML=J_mat;
Peaks = Peaks_1a;
[Peak_MSA,~,idy_MSA] = unique(Peaks,'rows');
C1 = tabulate(idy_MSA);
[B1,I1] = sort(C1(:,2),'descend');
Peak_MSA_1a = Peak_MSA(I1,:);
B1=B1/size(Peaks,1);
fre_1a =B1;
cum_fre_1a = cumsum(fre_1a);
Mutant_MSA_1a = zeros(length(B1),363);

Mutant_MSA = zeros(length(B1),length(ind_non_conserve));
for i = 1:length(B1)
    for j =1:length(ind_non_conserve)
        Start_site = num_mutants_combine_array_acc_all(j)+1;
        End_site = num_mutants_combine_array_acc_all(j+1);
        Mutant_MSA(i,j) = sum(Peak_MSA_1a(i,Start_site:End_site));
    end
end
Mutant_MSA_1a(:,ind_non_conserve) = Mutant_MSA;


all_p = zeros(length(B1),1);
for i = 1:length(B1)
    N=363;
    n = sum(Mutant_MSA_1a(i,:));
    q = sum(Mutant_MSA_1a(i,all_sites-383));
    j = length(all_sites);
    for k =q:min(j,n)  
        all_p(i) = all_p(i)+(nchoosek(j,k)*nchoosek(363-j,n-k))/nchoosek(N,n);
    end
end
all_p_1a = -log10(all_p);
all_p_1a(isinf(all_p_1a))=0;









load('Peaks.mat')
load('Model_1b.mat', 'num_mutants_combine_array');
load('Model_1b.mat', 'conserved')
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
Peaks = Peaks_1b;
[Peak_MSA,~,idy_MSA] = unique(Peaks,'rows');
C1 = tabulate(idy_MSA);
[B1,I1] = sort(C1(:,2),'descend');
Peak_MSA_1b = Peak_MSA(I1,:);
B1=B1/size(Peaks,1);
fre_1b =B1;
cum_fre_1b = cumsum(fre_1b);


ind_non_conserve = ones(1,363);
ind_non_conserve(conserved)=0;
ind_non_conserve = find(ind_non_conserve);
Mutant_MSA_1b = zeros(length(B1),363);

Mutant_MSA = zeros(length(B1),length(ind_non_conserve));
for i = 1:length(B1)
    for j =1:length(ind_non_conserve)
        Start_site = num_mutants_combine_array_acc_all(j)+1;
        End_site = num_mutants_combine_array_acc_all(j+1);
        Mutant_MSA(i,j) = sum(Peak_MSA_1b(i,Start_site:End_site));
    end
end
Mutant_MSA_1b(:,ind_non_conserve) = Mutant_MSA;


all_p = zeros(length(B1),1);
for i = 1:length(B1)
    N=363;
    n = sum(Mutant_MSA_1b(i,:));
    q = sum(Mutant_MSA_1b(i,all_sites-383));
    j = length(all_sites);
    for k =q:min(j,n)  
        all_p(i) = all_p(i)+(nchoosek(j,k)*nchoosek(363-j,n-k))/nchoosek(N,n);
    end
end
all_p_1b = -log10(all_p);
all_p_1b(isinf(all_p_1b))=0;






FIG=figure;


loglog(1:length(cum_fre_1a),cum_fre_1a,'Color',purple,'LineWidth',1);
hold on;


loglog(1:length(cum_fre_1b),cum_fre_1b,'Color',orange,'LineWidth',1);
yt = get(gca, 'YTick');
% plot([10 10], [min(yt) cum_fre_1b(10)], '--k','Color','k','LineWidth',1)
% plot([1 10], [cum_fre_1b(10) cum_fre_1b(10)], '--k','Color','k','LineWidth',0.7)
% plot([10 10], [min(yt) cum_fre_1a(10)], '--k','Color','k','LineWidth',0.7)
% plot([1 10], [cum_fre_1a(10) cum_fre_1a(10)], '--k','Color','k','LineWidth',0.7)
FIG.Name = 'log_cum_sum'
xlabel('Top x peaks');ylabel({'Cumulative fraction','of sequences in peak'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
set(gca,'TickLength',[0.02, 0.01])
FIG.Units = 'centimeters';
set(gca,'TickDir','out')
set(gcf,'Position',[10 10 4 5]);

set(gca,'Position',[.25 .25 .72 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.2 .2 .78 .74]);  %调整 XLABLE和YLABLE不会被切掉
% text(0,cum_fre_1a(10),char(vpa(cum_fre_1a(10),2)),'FontSize',8)
% text(0.34,cum_fre_1b(10),char(vpa(cum_fre_1b(10),2)),'FontSize',8)
set(gca,'YTick', [0:0.1:1])
yticklabels({'','0.1','','','','','','','','','1'})
set(gca,'XTick', [1 10 100])
xticklabels({'10^{0}','10^{1}','10^{2}'})
box off
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(gca,'TickLength',[0.035, 0.03])
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.25, 0.53, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 
% catch   end




FIG=figure;


hold on;

% l = categorical({'1a','1b'});
% l = reordercats(l,{'1a','1b'});
% b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% b.FaceColor = 'flat';
% b.CData(1,:)=purple;
% b.CData(2,:)=orange;
bar(1:2,[length(cum_fre_1a) zeros(1,1)], 0.2, 'FaceColor',purple,'LineWidth',0.2);
hold on
bar(1:2,[zeros(1,1) length(cum_fre_1b)], 0.2, 'FaceColor',orange,'LineWidth',0.2);
xlim([0.85 2.5])
set(gca,'XTick',1:2,'XTickLabel',...
    {'1a','1b'});
set(gca,'TickDir','out')
FIG.Name = 'barplot_num_peaks'
set(gca,'TickLength',[0.02, 0.01])
ylabel({'Number of peaks'})
% legend('Subtype 1a', 'Subtype 1b','Location','best')
FIG.Units = 'centimeters';


set(gcf,'Position',[5 5 3 5]);
% set(gcf,'Resize','off');
set(gca,'Position',[.78 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉

figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.32, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'TickLength',[0.035, 0.03])
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch     end
% 













FIG=figure;
enrich_1a = sum(Mutant_MSA_1a(:,all_sites-383),2)./sum(Mutant_MSA_1a,2);
enrich_1b = sum(Mutant_MSA_1b(:,all_sites-383),2)./sum(Mutant_MSA_1b,2);

% plot(enrich_1a,all_p_1a,'o','MarkerEdgeColor','None','MarkerFaceColor',purple,'MarkerSize',4)
scatter(enrich_1a,all_p_1a,'o','MarkerEdgeColor','None','MarkerFaceColor',purple,'SizeData',8);
hold on;
% plot(enrich_1b,all_p_1b,'o','MarkerEdgeColor','None','MarkerFaceColor',orange,'MarkerSize',5);
scatter(enrich_1b,all_p_1b,'o','MarkerEdgeColor','None','MarkerFaceColor',orange,'SizeData',8,'MarkerFaceAlpha',0.5);
percent_1a=sum(all_p_1a>-log10(0.05))/length(all_p_1a)
percent_1b=sum(all_p_1b>-log10(0.05))/length(all_p_1b)
plot(0:0.01:1, ones(size(0:0.01:1))*-log10(0.05), '--k','LineWidth',0.3)
xlim([0 1])
% set(gca,'YTick', [0 1 2 3])
% yticklabels({'0','1','2','3'})
set(gca,'TickDir','out')
set(gca,'TickLength',[0.035, 0.03])
% legend('Subtype 1a', 'Subtype 1b','Location','best')
% title('L2 reg. para. + 50%')
xlabel('Enrichment')
ylabel({'-log10(p-value)'})
% title(['Average epsilon=' num2str(single_error)])
FIG.Name = 'enrich_p_compare';
box off
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 4 5]);

set(gca,'Position',[.25 .25 .72 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.2 .2 .78 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.25, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch    end
% 

% %%
% FIG = figure;
% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',8)
% set(0,'DefaultTextFontSize',8)
% 
% 
% load('Peak_trajectory.mat')
% 
% percent_1a = sum(stat_1a)/sum(sum(stat_1a));
% percent_1a(2) = percent_1a(2)+percent_1a(3);
% percent_1a = percent_1a(1:2);
% percent_1a = percent_1a([2 1]);
% 
% percent_1b = sum(stat_1b)/sum(sum(stat_1b));
% percent_1b(2) = percent_1b(2)+percent_1b(3);
% percent_1b = percent_1b(1:2);
% percent_1b = percent_1b([2 1]);
% 
% b1=bar([percent_1a;zeros(1,2)]*100, 'stacked');
% for bars =b1
%     bars.BarWidth=0.2;
%     bars.FaceAlpha =1;
%     bars.LineWidth = 0.2;
% end
% 
% b1(1,1).FaceColor= purple;
% b1(1,2).FaceColor= purple;
% 
% 
% 
% b1(1,2).FaceAlpha= 0.4;
% 
% hold on 
% b2=bar([zeros(1,2);percent_1b]*100, 'stacked');
% for bars =b2
%     bars.BarWidth=0.2;
%     bars.FaceAlpha =1;
%      bars.LineWidth = 0.2;
% end
% 
% b2(1,1).FaceColor= orange;
% b2(1,2).FaceColor= orange;
% 
% b2(1,2).FaceAlpha= 0.4;
% 
% xtick =set(gca,'XTick',1:2,'XTickLabel',...
%     {'1a','1b'},'FontName','Arial');
% % h = legend('Self','Other top 10 peaks','Other peaks','Location','northoutside','Orientation','horizontal');
% % legend boxoff
% ylabel({'Percentage of trajectories, %'})
% box off
% xlim([0.85 2.5])
% set(gca,'XTick',1:2,'XTickLabel',...
%     {'1a','1b'});
% set(gca,'TickDir','out')
% FIG.Name = 'Peak_tract'
% set(gca,'TickLength',[0.02, 0.01])
% % ylabel({'Number of peaks'})
% % legend('Subtype 1a', 'Subtype 1b','Location','best')
% FIG.Units = 'centimeters';
% 
% 
% set(gcf,'Position',[5 5 3 5]);
% % set(gcf,'Resize','off');
% set(gca,'Position',[.78 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% % set(gcf,'Position',[10 10 7.84 6]);
% % set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'YTick', [0 50 100])
% figure_FontSize=8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(gca,'TickDir','out')
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.32, 0.5, 0]);
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% % set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% set(gca,'TickLength',[0.035, 0.03])
% % try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch        end
% % %%
% % 
% % set(0,'DefaultAxesFontName','Arial')
% % set(0,'DefaultTextFontName','Arial')
% % set(0,'DefaultAxesFontSize',8)
% % set(0,'DefaultTextFontSize',8)
% % 
% % 
% % load('Peak_trajectory.mat')
% % 
% % percent_1a = sum(stat_1a)/sum(sum(stat_1a));
% % percent_1a(2) = percent_1a(2)+percent_1a(3);
% % percent_1a = percent_1a(1:2);
% % 
% % 
% percent_1b = sum(stat_1b)/sum(sum(stat_1b));
% percent_1b(2) = percent_1b(2)+percent_1b(3);
% percent_1b = percent_1b(1:2);
% 
% 
% 
% FIG=figure;
% 
% 
% hold on;
% 
% % l = categorical({'1a','1b'});
% % l = reordercats(l,{'1a','1b'});
% % b = bar(l,[length(cum_fre_1a) length(cum_fre_1b)],0.2,'LineWidth',0.1,'BarWidth',1);
% % b.FaceColor = 'flat';
% % b.CData(1,:)=purple;
% % b.CData(2,:)=orange;
% bar(1:2,[percent_1a(1)*100 zeros(1,1)], 0.2, 'FaceColor',purple,'LineWidth',0.2);
% hold on
% bar(1:2,[zeros(1,1) percent_1b(1)*100], 0.2, 'FaceColor',orange,'LineWidth',0.2);
% xlim([0.85 2.5])
% set(gca,'XTick',1:2,'XTickLabel',...
%     {'1a','1b'});
% set(gca,'TickDir','out')
% FIG.Name = 'peak_tract'
% set(gca,'TickLength',[0.02, 0.01])
% ylabel({'Percentage of',' sequence trajectories, %'})
% % legend('Subtype 1a', 'Subtype 1b','Location','best')
% FIG.Units = 'centimeters';
% 
% 
% set(gcf,'Position',[5 5 3 5]);
% % set(gcf,'Resize','off');
% set(gca,'Position',[.78 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'Position',[.22 .25 .88 .71]);  %调整 XLABLE和YLABLE不会被切掉
% % set(gcf,'Position',[10 10 7.84 6]);
% % set(gca,'Position',[.18 .17 .76 .74]);  %调整 XLABLE和YLABLE不会被切掉
% 
% figure_FontSize=8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% % set(gca,'TickDir','out')
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.32, 0.5, 0]);
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% ylim([0 60])
% set(gca,'YTick', [0 30 60])
% set(gca,'TickLength',[0.035, 0.03])
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      end

