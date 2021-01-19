exposed_1a = [421,422,423,425,427,428,429,430,431,433,434,435,436,438,439,442,443,444,446,447,448,449,450,451,452,492,493,495,496,498,500,501,510,512,513,514,515,517,520,521,522,523,524,525,527,528,529,531,532,533,534,538,540,541,548,550,556,557,558,559,560,561,562,563,565,566,567,570,571,573,574,575,576,577,579,580,583,584,585,597,598,599,604,605,606,610,612,613,615,619,620,622,623,624,625,626,627,628,629,630,632,633,634,635,636,637,639,641,645];
%NC paper
exposed_1b=[410,411,412,414,415,416,417,418,419,421,422,423,424,427,428,430,431,432,433,434,435,438,439,442,443,444,445,446,447,448,449,450,451,453,454,455,456,457,458,460,461,463,464,466,467,469,470,471,473,474,476,477,478,479,480,481,482,483,484,485,488,489,490,492,493,496,498,501,510,512,513,514,515,520,521,522,523,524,525,526,527,528,529,531,532,533,534,538,540,542,543,545,546,548,549,556,558,560,561,562,570,573,574,576,577,578,584,586,587,590,591,593,595,596,604,606,610,612,622,623,624,625,626,627,628,629,630,632,634,635,636,637,639,641,645,646,647];
[~,ia,ib] = intersect(exposed_1a,exposed_1b);

exposed_1b = exposed_1b(ib);
exposed_1a = exposed_1a(ia);

load('escape_time_1a_500.mat', 'mean_escape_time');
load('mean_escape_time_1b.mat')
escape_time_1b = all_mean_escape_time_1b(exposed_1b-383);
escape_time_1a = mean_escape_time(exposed_1a-383);

% label_xaxis_data = cellstr(num2str(exposed_1a'));

label_xaxis_data = cellstr(num2str(exposed_1a'));


% show the first half
label_xaxis_data  = label_xaxis_data (1:40);
escape_time_1a = escape_time_1a(1:40);
escape_time_1b= escape_time_1b(1:40);
thresh=200;
idx_high = escape_time_1a>=thresh & escape_time_1b>=thresh;

label_yaxis_data = {'Subtype 1a','Subtype 1b'};



% label_yaxis_data = {80,100,120,140};



map = zeros(length(label_yaxis_data),length(label_xaxis_data));
for kk = 1:length(escape_time_1a)
    map(1,kk) = floor(escape_time_1a(kk));
end
for kk = 1:length(escape_time_1b)
    map(2,kk) = floor(escape_time_1b(kk));
end

% map  = flip(map,1);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

FIG=figure;
% subplot(2,1,2)
% [grad,im]=colorGradient(	hex2rgb('#FCFAF2'),hex2rgb('#4A225D'),128);
% [grad,im]=colorGradient(	[0.21569,0.49412,0.72157],[0.89412,0.10196,0.1098],128);
load('color_gradient.mat')
grad = gradient;
h=heatmap(label_xaxis_data,label_yaxis_data,map,'GridVisible','off','Colormap',grad,'ColorLimits',[0 200],'ColorbarVisible','off','FontName','Arial','FontSize',8);

% set(struct(h).NodeChildren(3), 'XTickLabelRotation', 90);
% text('Units', 'Normalized', 'Position', [-0.09,0.5],'String','Escape time threshold, \tau','Rotation',90)
% ,'FontName','Arial','FontSize',8,'FontWeight','Bold'

FIG.Units = 'centimeters';
FIG.Name = 'RBall_thre_1b';
set(gcf,'Position',[10 10 20 2.4]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.08 .3 .91 .6]);  %调整 XLABLE和YLABLE不会被切掉

% set(FIG, 'DefaultTextFontSize', 8);
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top','FontWeight','Bold');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle','FontWeight','Bold');
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',0);
% set(gca, 'LineWidth',1)
% set(gca,'FontName','Arial','FontSize',8,'FontWeight','Bold')


% text(0.2,0.6,'Escape time threshold, \tau','Rotation',90)
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12  12])
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end

FIG=figure

hold on

idx_high = find(idx_high);

for i = idx_high
label_xaxis_data{i,1} = (convertStringsToChars(strjoin(['\bf',label_xaxis_data(i)])));
text(i-0.2,2,'\ast','FontSize',12)
end

FIG.Name  = 'xlabel';

set(gca,'XTick',1:40,'XTickLabel',label_xaxis_data);

set(gca,'YTick',[])
xlim([0.5 40.5])
FIG.Units = 'centimeters';
FIG.Name = 'xlabel';
set(gcf,'Position',[10 10 20 2.4]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.08 .3 .91 .3]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'YTick',[]);

xtickangle(45)
% ylabel('Minimum escape time, t^{anti}')
% ylabel('Minimum escape time, t^{min}_{e}')


% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end

% print(FIG, '-depsc2', '-r600', '-tiff', '-loose', ['C:\Users\27909\Desktop\' FIG.Name]);

%% show the second half
exposed_1a = [421,422,423,425,427,428,429,430,431,433,434,435,436,438,439,442,443,444,446,447,448,449,450,451,452,492,493,495,496,498,500,501,510,512,513,514,515,517,520,521,522,523,524,525,527,528,529,531,532,533,534,538,540,541,548,550,556,557,558,559,560,561,562,563,565,566,567,570,571,573,574,575,576,577,579,580,583,584,585,597,598,599,604,605,606,610,612,613,615,619,620,622,623,624,625,626,627,628,629,630,632,633,634,635,636,637,639,641,645];
%NC paper
exposed_1b=[410,411,412,414,415,416,417,418,419,421,422,423,424,427,428,430,431,432,433,434,435,438,439,442,443,444,445,446,447,448,449,450,451,453,454,455,456,457,458,460,461,463,464,466,467,469,470,471,473,474,476,477,478,479,480,481,482,483,484,485,488,489,490,492,493,496,498,501,510,512,513,514,515,520,521,522,523,524,525,526,527,528,529,531,532,533,534,538,540,542,543,545,546,548,549,556,558,560,561,562,570,573,574,576,577,578,584,586,587,590,591,593,595,596,604,606,610,612,622,623,624,625,626,627,628,629,630,632,634,635,636,637,639,641,645,646,647];
[~,ia,ib] = intersect(exposed_1a,exposed_1b);

exposed_1b = exposed_1b(ib);
exposed_1a = exposed_1a(ia);

load('escape_time_1a_500.mat', 'mean_escape_time');
load('mean_escape_time_1b.mat')
escape_time_1b = all_mean_escape_time_1b(exposed_1b-383);
escape_time_1a = mean_escape_time(exposed_1a-383);

% label_xaxis_data = cellstr(num2str(exposed_1a'));

label_xaxis_data = cellstr(num2str(exposed_1a'));
label_xaxis_data  = label_xaxis_data (41:end);
escape_time_1a = escape_time_1a(41:end);
escape_time_1b= escape_time_1b(41:end);
thresh=200;
idx_high = escape_time_1a>=thresh & escape_time_1b>=thresh;

label_yaxis_data = {'Subtype 1a','Subtype 1b'};



% label_yaxis_data = {80,100,120,140};



map = zeros(length(label_yaxis_data),length(label_xaxis_data));
for kk = 1:length(escape_time_1a)
    map(1,kk) = floor(escape_time_1a(kk));
end
for kk = 1:length(escape_time_1b)
    map(2,kk) = floor(escape_time_1b(kk));
end

% map  = flip(map,1);
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

FIG=figure;
% subplot(2,1,2)
% [grad,im]=colorGradient(	hex2rgb('#FCFAF2'),hex2rgb('#4A225D'),128);
% [grad,im]=colorGradient(	[0.21569,0.49412,0.72157],[0.89412,0.10196,0.1098],128);
load('color_gradient.mat')
grad = gradient;
h=heatmap(label_xaxis_data,label_yaxis_data,map,'GridVisible','off','Colormap',grad,'ColorLimits',[0 200],'ColorbarVisible','off','FontName','Arial','FontSize',8);

% set(struct(h).NodeChildren(3), 'XTickLabelRotation', 90);
% text('Units', 'Normalized', 'Position', [-0.09,0.5],'String','Escape time threshold, \tau','Rotation',90)
% ,'FontName','Arial','FontSize',8,'FontWeight','Bold'

FIG.Units = 'centimeters';
FIG.Name = 'RBall_thre_1b';
set(gcf,'Position',[10 10 20 2.4]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.08 .3 .91 .6]);  %调整 XLABLE和YLABLE不会被切掉

% set(FIG, 'DefaultTextFontSize', 8);
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top','FontWeight','Bold');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle','FontWeight','Bold');
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',0);
% set(gca, 'LineWidth',1)
% set(gca,'FontName','Arial','FontSize',8,'FontWeight','Bold')


% text(0.2,0.6,'Escape time threshold, \tau','Rotation',90)
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 12  12])
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end

FIG=figure

hold on

idx_high = find(idx_high);

for i = idx_high
label_xaxis_data{i,1} = (convertStringsToChars(strjoin(['\bf',label_xaxis_data(i)])));
text(i-0.2,2,'\ast','FontSize',12)
end

FIG.Name  = 'xlabel';

set(gca,'XTick',1:40,'XTickLabel',label_xaxis_data);

set(gca,'YTick',[])
xlim([0.5 40.5])
FIG.Units = 'centimeters';
FIG.Name = 'xlabel';
set(gcf,'Position',[10 10 20 2.4]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.08 .3 .91 .3]);  %调整 XLABLE和YLABLE不会被切掉
% set(gca,'YTick',[]);

xtickangle(45)
% ylabel('Minimum escape time, t^{anti}')
% ylabel('Minimum escape time, t^{min}_{e}')


try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end

% print(FIG, '-depsc2', '-r600', '-tiff', '-loose', ['C:\Users\27909\Desktop\' FIG.Name]);

%% figure 5b can be generated by the respective .pse file