run startup.m
%1a
load('Model_1a.mat', 'msa_aa')
load('Model_1a.mat', 'weight_id')
load('Model_1a.mat', 'ind_non_conserve')
sequences = mat2cell(msa_aa,ones(1,size(msa_aa,1)),size(msa_aa,2));
[Profile, Symbols] = seqprofile(sequences);
site_entropy_1a=zeros(1,363);
entropy = Profile.*log(Profile);
for i=1:363
    

    site_entropy_1a(i) = -sum(entropy(~isnan(entropy(:,i)),i));
end

% 1b
load('Model_1b.mat', 'sequences')
load('Model_1b.mat', 'weight_seq')
load('Model_1b.mat', 'conserved')
[Profile, Symbols] = seqprofile(sequences);
site_entropy_1b=zeros(1,363);
entropy = Profile.*log(Profile);
for i=1:363
    

    site_entropy_1b(i) = -sum(entropy(~isnan(entropy(:,i)),i));
end




data = [site_entropy_1a site_entropy_1b];
G = [zeros(size(site_entropy_1a)) ones(size(site_entropy_1b))];


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
median_lineWidth = 1.5;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = '';
outlier_markerSize = 3.5;
outlier_marker_edgeWidth = 0.001;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'Subtype 1a',sprintf('Subtype 1b')};
text_ylabel = 'Entropy at each residue';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 5;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')
x1=0.8+0.4*(rand(length(site_entropy_1a),1));
x2=1.8+0.4*(rand(length(site_entropy_1b),1));
size_marker = 20;
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
f1=scatter(x1,site_entropy_1a ,'o','MarkerEdgeColor','w','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,site_entropy_1b,'o','MarkerEdgeColor','w','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on




figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);
% V=violinplot(data, [repmat("Subtype 1a",1,length(site_entropy_1a)) repmat("Subtype 1b",1,length(site_entropy_1b))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel('Entropy of each residue');
% V(1, 1).ViolinColor = purple;
% V(1, 2).ViolinColor = orange;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;

set(gca,'YTick',0:1:2)
yt = get(gca, 'YTick');
axis([xlim    0  2.5])
xt = get(gca, 'XTick');
hold on
% plot(xt([1 2]), [1 1]*0.25, '-k','LineWidth',0.3)
% plot(xt([1 1]), [0.95 1]*0.25, '-k','LineWidth',0.3)
% plot(xt([2 2]), [0.95 1]*0.25, '-k','LineWidth',0.3)
% P = ranksum(mutantperseq_1b,mutantperseq_1a)
% text(0.9,0.27,'$$ P = 1.69 \times 10^{-31} $$','interpreter','latex','FontSize',12)
yt = get(gca, 'YTick');
% axis([xlim    0  ceil(max(yt)*1.2)])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*2.5, '-k','LineWidth',0.5)
% plot(xt([1 1]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
% plot(xt([2 2]), [0.95 1]*max(yt)*1.05, '-k','LineWidth',0.5)
P = ranksum(site_entropy_1b,site_entropy_1a,'tail','right')
% text(1.15,max(yt)*1.1,'$$ P = 0.0928 $$','interpreter','latex','FontSize',14)
% text(0.8,-100,'residues','FontSize',8)



FIG.Name = 'Mutant_seq_compare';


FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 10 8]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(gca,'TickDir','out')
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.01])
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
