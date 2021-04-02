%%1a
load('epitopes_1a.mat')
ind=[];

for i = 1:length(starting)
    for j=i+1:length(starting)
        common = length(intersect(starting(i):ending(i),starting(j):ending(j)));
        
        if common == length(starting(i):ending(i)) 
            ind = [ind; i];
        else
            if common == length(starting(j):ending(j))
                
            ind =[ind ;j];
            end
        end
        
    end
end
starting(ind)=[];
ending(ind)=[];
ind_con(ind)=[];
cov_1a= [];
for i = 1:length(starting)
    cov_1a = [cov_1a starting(i):ending(i)];
end

ind = zeros(1,length(pos));



for i =1:length(pos)
    if all(pos{1,i}>=384 & pos{1,i}<=746)
        ind(i)=1;
    end
    
end
ind = logical(ind);
peptides_dis = peptides_dis(ind);
pos = pos(ind);
ind_dis = ind_dis(ind);
ind=[];
for i = 1:length(pos)
    for j=i+1:length(pos)
        common = length(intersect(pos{1,i},pos{1,j}));
        
        if common == length(pos{1,i}) 
            ind = [ind; i];
        else
            if common == length(pos{1,j})
                
            ind =[ind ;j];
            end
        end
        
    end
end
pos(ind)=[];
ind_dis(ind)=[];
for i = 1:length(pos)
    cov_1a = [cov_1a pos{1,i}];
end
cov_1a = unique(cov_1a);

unique_1a = length(starting)+length(pos)
ind_1a = [ind_con ind_dis];
%% 1b
load('epitopes_1b.mat')
ind=[];
for i =1:size(ind_dis,1)
    tmp{i} = ind_dis(i,:);
end
ind_dis=tmp;


for i = 1:length(starting)
    for j=i+1:length(starting)
        common = length(intersect(starting(i):ending(i),starting(j):ending(j)));
        
        if common == length(starting(i):ending(i)) 
            ind = [ind; i];
        else
            if common == length(starting(j):ending(j))
                
            ind =[ind ;j];
            end
        end
        
    end
end
starting(ind)=[];
ending(ind)=[];
ind_con(ind)=[];
cov_1b= [];
for i = 1:length(starting)
    cov_1b = [cov_1b starting(i):ending(i)];
end

ind = zeros(1,length(pos));



for i =1:length(pos)
    if all(pos{1,i}>=384 & pos{1,i}<=746)
        ind(i)=1;
    end
    
end
ind = logical(ind);
peptides_dis = peptides_dis(ind);
pos = pos(ind);
ind_dis = ind_dis(ind);
ind=[];
for i = 1:length(pos)
    for j=i+1:length(pos)
        common = length(intersect(pos{1,i},pos{1,j}));
        
        if common == length(pos{1,i}) 
            ind = [ind; i];
        else
            if common == length(pos{1,j})
                
            ind =[ind ;j];
            end
        end
        
    end
end
pos(ind)=[];
ind_dis(ind)=[];
for i = 1:length(pos)
    cov_1b = [cov_1b pos{1,i}];
end
cov_1b = unique(cov_1b);

unique_1b = length(starting)+length(pos)
ind_1b = [ind_con ind_dis];

%% load ref_B_cell

load('IEDB_id.mat')
[~,I,miss]=intersect(cell2mat(epitope_id),str2double(ind_1a'));

%6 missing ones to check manually
% all from ref id: 1001984 
miss = str2double(ind_1a(setdiff(1:71,miss)));


ref_1a = unique(cell2mat(ref_id(I)));


[~,I,~]=intersect(cell2mat(epitope_id),str2double(ind_1b'));

ref_1b = unique(cell2mat(ref_id(I)));


%% plot figure
label_xaxis_data = cellstr(num2str((1:363)'));


label_yaxis_data = {'Subtype 1a','Subtype 1b'};
map = zeros(length(label_yaxis_data),length(label_xaxis_data));

for kk = 1:length(cov_1a)
    map(1,cov_1a(kk)-383) = 1;
end
for kk = 1:length(cov_1b)
    map(2,cov_1b(kk)-383) = 1;
end

FIG=figure;

heatmap(label_xaxis_data,label_yaxis_data,map,'GridVisible','off','Colormap',[1 1 1; [168 216 185]/255],'ColorLimits',[0 max(max(map))],'ColorbarVisible','off','FontName','Arial','FontSize',8);
FIG.Units = 'centimeters';
FIG.Name = 'epitopes';
set(gcf,'Position',[10 10 20 2.4]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.09 .3 .88 .6]);  %调整 XLABLE和YLABLE不会被切掉
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));

% text('Units', 'Normalized', 'Position', [0.08,0.2],'String','384')
% text('Units', 'Normalized', 'Position', [0.92,0.2],'String','384')
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
