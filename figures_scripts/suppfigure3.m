
load('DCA_Ahmed.mat');
num=363;
F = DI_new;
F =triu(F,1);
Nmax = num*(num-1)/2; % get Nmax biggest entries
[ ~, Ind ] = sort(F(:),1,'descend');
[ ind_row, ind_col ] = ind2sub(size(F),Ind(1:Nmax)); % fetch indices
load('6mei.mat')
indx=indx-383;
indy=indy-383;
ia1 = ind_row<=(647-383) & ind_row>= (410-383);
ia2 = ind_col<=(647-383) & ind_col>= (410-383);
ia = ia1 & ia2;
ind_col =ind_col(ia);
ind_row =ind_row(ia);


% load('I-TASSER\model5.mat')
TPR = zeros(size(ind_row));
count =0;
for i =1:length(ind_row)
   ia1 = find(indx==ind_row(i));
   ia2 = find(indy==ind_col(i));
   if(intersect(ia1,ia2))
        count = count+1;
   end

    TPR(i) = count/i;

end

FIG=figure
color = [0.6000    0.6000    0.6000];%gray
plot(TPR,'Color',color,'LineWidth',1);

% 
% box off
set(gca, 'XScale', 'log')
xlim([1 1e3])

xlabel('Top x pairs')
ylabel({'Precision'})
FIG.Units = 'centimeters';
FIG.Name = 'DCA';

FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 10 8]);
set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
box off
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'TickLength',[0.02, 0.03])
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
