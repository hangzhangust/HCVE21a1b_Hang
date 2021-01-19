function bh = figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm)

% Code for generating boxplot figures
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%% Main box plot
% figure;
bh = boxplot(data,G,...
    'whisker',whisker_value,'symbol',outlier_marker,...
    'color','k','jitter',outlier_jitter_value,...
    'labels',label_xaxis_data,...
    'widths',box_widths_value,'LabelOrientation',label_orientation_choice);

set(bh,'linewidth',box_lineWidth);
ylabel(text_ylabel)
title(text_title)

% if length(label_xaxis_data)>10 %rotate xaxis labels by 45 degrees and make their font smaller only if labels>10
%     set(gca,'XTickLabelRotation',45)
%     xL = xlabel(text_xlabel);
%     ax = ancestor(gca, 'axes');
%     xrule = ax.XAxis;
%     xrule.FontSize = 10;
%     xL.FontSize = 10;
% end


% Coloring each box
h = findobj(gca,'Tag','Box');
if size(box_color,1) ~= 1   %if colors provided for each box
    for kk = 1:length(h)
        patch(get(h(kk),'XData'),get(h(kk),'YData'),box_color(length(h)-kk+1,:),'FaceAlpha',box_color_transparency);
        %     h(kk).Color = 'w';color_boxplot(length(HmAb_Pierce2016)-kk+1,:);
    end
else
    for kk = 1:length(h)
        patch(get(h(kk),'XData'),get(h(kk),'YData'),box_color,'FaceAlpha',box_color_transparency);
    end
end

% Adjusting median
h=findobj(gca,'tag','Median');
for kk = 1:length(h)
    h(kk).LineWidth = median_lineWidth;
    h(kk).Color = median_color;
end

% Adjusting outliers
h=findobj(gca,'tag','Outliers');
for kk = 1:length(h)
    if size(box_color,1) ~= 1   %if colors provided for each box
        h(kk).MarkerFaceColor = box_color(length(h)-kk+1,:);
    else
        h(kk).MarkerFaceColor = box_color;
    end
    h(kk).MarkerEdgeColor = outlier_marker_edgeColor;
    h(kk).MarkerSize = outlier_markerSize;
    h(kk).LineWidth = outlier_marker_edgeWidth;
end

% alpha(h,0.7)

ylim([ylim_min ylim_max])
post_proc_fig_box

% Saving figure
if savefig == 1
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 fig_width_cm fig_height_cm])
    print(savefig_name,'-dpng','-r600')
end