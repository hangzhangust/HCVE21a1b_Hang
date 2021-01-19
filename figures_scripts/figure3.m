%% predcit known escape mutations using escape times
run startup.m
rng default
L=363;
load('mean_escape_time_1b.mat', 'all_mean_escape_time_1b')
all_E=all_mean_escape_time_1b;
% polymorphisms_associated_with_neutralization_resistance=[415;416;417;422;424;431;433;435;438;442;444;446;447;453;456;458;461;466;475;478;482;501;520;524;528;531;533;538;557;558;560;580;610;636;655;665;713];
% exclude 100% conserved site: 520 
polymorphisms_associated_with_neutralization_resistance=[384;386;388;390;391;393;394;...
    395;396;397;398;399;400;401;402;403;404;405;407;408;410;415;416;417;422;424;431;433;...
    434;435;438;442;444;446;453;456;461;466;475;482;501;524;528;531;533;538;557;558;560;580;608;610;636;713;];
%boxplot

G = [zeros(1,length(polymorphisms_associated_with_neutralization_resistance)) ...
    ones(1,length(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)))];
data = [all_E(polymorphisms_associated_with_neutralization_resistance-383) ...
    all_E(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383))];

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
L=363;
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'Escape' ,sprintf('Remaining')};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 500;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')
size_marker = 20;
x1=0.8+0.4*(rand(length(polymorphisms_associated_with_neutralization_resistance),1));
x2=1.8+0.4*(rand(length(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),1));
f1=scatter(x1,all_E(polymorphisms_associated_with_neutralization_resistance-383) ,'o','MarkerEdgeColor','w','MarkerFaceColor',blue,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,all_E(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),'o','MarkerEdgeColor','w','MarkerFaceColor',red,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on

hold on
figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);
hold on;


red = color_scheme_npg(1,:);
blue = color_scheme_npg(2,:);
% V=violinplot(data, [repmat("Escape",1,length(polymorphisms_associated_with_neutralization_resistance)) repmat("Remaining",1,length(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel('Escape time');
% V(1, 1).ViolinColor = blue;
% V(1, 2).ViolinColor = red;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;
% V(1, 1).BoxColor  = 'None';
% V(1, 1).BoxColor  = 'None';
set(gca,'YTick',0:225:450)
yt = get(gca, 'YTick');
axis([xlim    0  540])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*472.5, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
P = ranksum(all_E(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),all_E(polymorphisms_associated_with_neutralization_resistance-383),'tail','right');
ind = floor(log10(P));
P = roundn(P/10^ind,-1);
t = ['$$ P = ',num2str(P),' \times 10^{',num2str(ind),'} $$'];

text(1,500,t,'interpreter','latex','FontSize',8)
% text(1,max(yt)*1.15,'$$ P = 9.3 \times 10^{-20} $$','interpreter','latex','FontSize',8)
% text(0.8,-100,'residues','FontSize',8,'FontName', 'Arial')

set(gca,'TickLength',[0.02, 0.01])
FIG.Name = 'Escape_mutation';
set(gca,'TickDir','out')
FIG.Units = 'centimeters';
FIG.Name = 'Escape_compare';
set(gcf,'Position',[6.53 6.53 5 10]);
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'Position',[.25 .1 .72 .88]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.235, 0.5, 0]);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end

% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',1);
% set(gca,'FontName','Arial','FontSize',8)


% % export_fig C:\Users\27909\Desktop\Escape_mutation.png -native
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');


hold off

% fprintf('\nP = %.1e, Mann-Whitney test\n',P)

%% predict buried and exposed residues
color_scheme_lancet = [...
         0    0.2745    0.5451; ...
    0.9294         0         0; ...
    0.2588    0.7098    0.2510; ...
         0    0.6000    0.7059; ...
    0.5725    0.3686    0.6235; ...
    0.9922    0.6863    0.5686; ...
    0.6784         0    0.1647; ...
    0.6784    0.7137    0.7137; ...
    0.1059    0.0980    0.0980];
run startup.m
L=363;
load('mean_escape_time_1b.mat', 'all_mean_escape_time_1b')
all_E=all_mean_escape_time_1b;
exposed_residues=[410,411,412,414,415,416,417,418,419,421,422,423,424,427,428,430,431,432,433,434,435,438,439,442,443,444,445,446,447,448,449,450,451,453,454,455,456,457,458,460,461,463,464,466,467,469,470,471,473,474,476,477,478,479,480,481,482,483,484,485,488,489,490,492,493,496,498,501,510,512,513,514,515,520,521,522,523,524,525,526,527,528,529,531,532,533,534,538,540,542,543,545,546,548,549,556,558,560,561,562,570,573,574,576,577,578,584,586,587,590,591,593,595,596,604,606,610,612,622,623,624,625,626,627,628,629,630,632,634,635,636,637,639,641,645,646,647];
buried_residues = setdiff(410:647,exposed_residues);
Buried=all_E(buried_residues-383);
% Exposed= all_E(~Bury);
Exposed= all_E(exposed_residues-383);

group = [zeros(size(Exposed)),ones(size(Buried))];
% boxplot([B,Ex],group,'Labels',{'Buried','Exposed'});
green = color_scheme_npg(3,:);
lightgreen = color_scheme_lancet(7,:);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
box_lineWidth = 0.75;
box_widths_value = 0.5;
black = [0 0 0];
box_color = [black;black];
box_color_transparency = 0; %faceAlpha
median_lineWidth = 0.75;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;label_xaxis_data = {'Exposed','Buried'};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 500;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;

FIG=figure;
x1=0.8+0.4*(rand(length(Exposed),1));
x2=1.8+0.4*(rand(length(Buried),1));
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
f1=scatter(x1,Exposed ,'o','MarkerEdgeColor','w','MarkerFaceColor',lightgreen,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,Buried,'o','MarkerEdgeColor','w','MarkerFaceColor',green,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on
% f1=scatter(x1,Exposed ,'o','MarkerEdgeColor','none','MarkerFaceColor',lightgreen,'SizeData',size_marker*0.8,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
% f1=scatter(x1,Exposed ,'o','MarkerEdgeColor','w','MarkerFaceColor',lightgreen,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 



hold on;
figure_boxplot([Exposed,Buried],group,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);
hold on

% V=violinplot([Exposed,Buried], [repmat("Exposed",1,length(Exposed)) repmat("Buried",1,length(Buried))],'GroupOrder',{'Exposed','Buried'},'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel('Escape time');
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
% 
% V(1, 1).ViolinColor = lightgreen;
% V(1, 2).ViolinColor = green;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;
set(gca,'YTick',0:225:450)
yt = get(gca, 'YTick');
axis([xlim    0  540])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*472.5, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*525, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*525, '-k','LineWidth',0.5)
P = ranksum(Buried,Exposed,'tail','right');

ind = floor(log10(P));
P = roundn(P/10^ind,-1);
t = ['$$ P = ',num2str(P),' \times 10^{',num2str(ind),'} $$'];

text(1,500,t,'interpreter','latex','FontSize',8)
FIG.Units = 'centimeters';
FIG.Name = 'Bury_Exposed';
set(gca,'TickDir','out')
set(gcf,'Position',[6.53 6.53 5 10]);
set(gca,'TickLength',[0.02, 0.01])
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'Position',[.25 .1 .72 .88]);  %调整 XLABLE和YLABLE不会被切掉
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.235, 0.5, 0]);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end

