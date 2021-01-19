 %%
 load('protein_entropy.mat')

 FIG=figure

 set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
scale = 3;
names = {'Core';'E1';'\color[rgb]{0.7961 0.1059 0.2706}E2';'p7';'NS2';'NS3';'NS4A';'NS4B';'NS5A';'NS5B'}; %#cb1b45
scatter( all_len, -all_entropy, all_num/10^4*scale,[0.6 0.6 0.6],'filled');

hold on;


scatter( 363, -entropy_1a, num_of_para_1a/10^4*scale,[0.6 0.6 0.6],'filled');

text( 363+1, -entropy_1a-0.015, '\color[rgb]{ 0.6235  0.2078 0.2275}E2 (1a)') %#9f351c

a = gca;
% set box property to off and remove background color
set(a,'box','off','color','none')
% create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
% set original axes as active
axes(a)
%     % link axes in case of zooming
%     linkaxes([a b])
 dx = ones(10,1);
 dy = ones(10,1)*-0.015;
%  dy = dy-[0.01;0.03;0.04];
%  dx = dx-[0;2;2;];
 
 
 
 text( all_len+dx, -all_entropy+dy, names)
 xlabel('Protein length')
 ylabel('Mean residue entropy')
%  
% text( 600, 0.25, 'Num of parameters')
hold on;  h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ok','MarkerFaceColor','black','MarkerSize',sqrt(scale));
h(2) = plot(NaN,NaN,'ok','MarkerFaceColor','black','MarkerSize',sqrt(scale*10));
h(3) = plot(NaN,NaN,'ok', 'MarkerFaceColor','black','MarkerSize',sqrt(scale*100));
ld =legend(h, '1e+04','1e+05','1e+06');
legend boxoff;

pos = [0.71 0.597142857142858 0.145297619626636 0.297460356333903];
annotation(FIG,'rectangle',...
    pos);
set(ld,...
    'Position',pos,...
    'FontSize',8);
title(ld,['Number of',sprintf('\n'),'parameters']);

set(gca,'TickDir','out'); 
set(get(gca,'XLabel'),'FontSize',8,'Vertical','top');
% set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',8,'Vertical','middle');
grid on
FIG.Name = 'protein_length'
FIG.Units = 'centimeters';
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);

% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');