%% compare escape time
load('mean_escape_time_1b.mat')
load('escape_time_1a_500.mat')
run startup.m

rng default

exposed_1a = [421,422,423,425,427,428,429,430,431,433,434,435,436,438,439,442,443,444,446,447,448,449,450,451,452,492,493,495,496,498,500,501,510,512,513,514,515,517,520,521,522,523,524,525,527,528,529,531,532,533,534,538,540,541,548,550,556,557,558,559,560,561,562,563,565,566,567,570,571,573,574,575,576,577,579,580,583,584,585,597,598,599,604,605,606,610,612,613,615,619,620,622,623,624,625,626,627,628,629,630,632,633,634,635,636,637,639,641,645];
%NC paper
exposed_1b=[410,411,412,414,415,416,417,418,419,421,422,423,424,427,428,430,431,432,433,434,435,438,439,442,443,444,445,446,447,448,449,450,451,453,454,455,456,457,458,460,461,463,464,466,467,469,470,471,473,474,476,477,478,479,480,481,482,483,484,485,488,489,490,492,493,496,498,501,510,512,513,514,515,520,521,522,523,524,525,526,527,528,529,531,532,533,534,538,540,542,543,545,546,548,549,556,558,560,561,562,570,573,574,576,577,578,584,586,587,590,591,593,595,596,604,606,610,612,622,623,624,625,626,627,628,629,630,632,634,635,636,637,639,641,645,646,647];
% 
% [buried_residues,exposed_1b]=classify_buried_exposed('raptorx.txt');
mean_escape_time = mean_escape_time(exposed_1a-383);
all_mean_escape_time_1b = all_mean_escape_time_1b(exposed_1b-383);


data = [mean_escape_time  all_mean_escape_time_1b ];
G = [zeros(size(mean_escape_time)) ones(size(all_mean_escape_time_1b))];


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
label_xaxis_data = {'1a',sprintf('1b')};
text_ylabel = {'Escape time of',' exposed residues'};
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
outlier_jitter_value = 0;
savefig = 0;
savefig_name = 'escape_time_compare';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')
x1=0.8+0.4*(rand(length(mean_escape_time),1));
x2=1.8+0.4*(rand(length(all_mean_escape_time_1b),1));
size_marker = 10;
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
f1=scatter(x1,mean_escape_time ,'o','MarkerEdgeColor','k','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,all_mean_escape_time_1b,'o','MarkerEdgeColor','k','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on



figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

% V=violinplot(data, [repmat("1a",1,length(mean_escape_time)) repmat("1b",1,length(all_mean_escape_time_1b))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel({'Escape time of',' exposed residues'});
% V(1, 1).ViolinColor = purple;
% V(1, 2).ViolinColor = orange;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;


set(gca,'YTick',0:250:500)
yt = get(gca, 'YTick');
axis([xlim    0  600])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*550, '-k','LineWidth',0.3)
% plot(xt([1 1]), [0.95 1]*550, '-k','LineWidth',0.3)
% plot(xt([2 2]), [0.95 1]*550, '-k','LineWidth',0.3)
P = ranksum(all_mean_escape_time_1b,mean_escape_time,'tail','left')

ind = floor(log10(P));
P = roundn(P/10^ind,-1);
t = ['$$ P = ',num2str(P),' \times 10^{',num2str(ind),'} $$'];

% text(0.9,600,t,'interpreter','latex','FontSize',8,'FontName','Arial')


FIG.Name = 'escape_time_compare';
set(gca,'fontname','Arial')  % Set it to times

FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 8 8.58]);
set(gcf,'Position',[6.53 6.53 4 4.16]);
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'TickLength',[0.035, 0.01])
set(gca,'Position',[.37 .12 .6 .86]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.01])
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.5, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% 
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end

%% figure 4 b can be generated by the respective .pse file

%% Box plot of minumum escape times of binding sites of each antibody
load('escape_time_1a_500.mat', 'mean_escape_time')
all_E=mean_escape_time;
HmAbs_all=[];

HmAb_Pierce2016 = [];

%CBH-4D
HmAb_Pierce2016{1} = [494 497 502 504 506:509 511 537 539 540 542:545 547 549:552 554 556 559 561 562 564 565 584 585 592 594 598 600 602 603 607:611 614 617:619 621 623 624 626 627 629:633 638 640 642:644];
%CBH-4G
HmAb_Pierce2016{2} = [494 497 503 505:509 511 517 537 539 540 542:545 547 549:552 554 556 559 561 564 565 584 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 629 631:633 638 640 642:644];
%CBH-4B
HmAb_Pierce2016{3} = [494 497 503 505:509 537 539 540 550:552 554 559 561 564 565 600 602 603 607:611 614 617:619 621 624 627 629 631:633 638 640 642:644];
%CBH-20
% HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 611 614 617:619 621 623 624 627 631 632 638 640 642:644];
%included the residues with close to threshold (RB = 21,22) as well
HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 610 611 614 617:619 621 623 624 627 631 632 638 640 642:644]; 
%CBH-21
HmAb_Pierce2016{5} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];
%CBH-22
HmAb_Pierce2016{6} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];

%HC-1
HmAb_Pierce2016{7} = [429 494 503:506 508 509 529 530 535 537 539 552 554 559 564 607 611 614 617 644];
%HC-11
HmAb_Pierce2016{8} = [425 428 429 436:438 442 443 494 497 502:504 506:509 511 520 530 535 537 539 550:552 554 556 559 564 565 602 603 607 608 611 614 617:621 624 640 643 644];
%A27
HmAb_Pierce2016{9} = [424:429 437 438 494 497 499 502:504 506:509 511 520 529 530 535 537 539 540 550:552 554 556 559 564 565 602 603 607:611 614 616:619 621 623 624 638 640 642:644];

%CBH-23
HmAb_Pierce2016{10} = [494 508 509 537 539 549 552 554 564 611 614 644]; %same as CBH7
%CBH-7
HmAb_Pierce2016{11} = [494 506 508 509 537 539 549 552 554 564 611 614 617 621 644];

%HC84-20
HmAb_Pierce2016{12} = [429 441 494 497 502:509 511 537 539 552 554 559 564 607 608 611 613 614 616:619 621 640 643 644];
%HC84-24
HmAb_Pierce2016{13} = [429 442 443 494 497 502:509 511 537 539 552 554 559 564 607 608 611 614 617:619 621 643 644];
%HC84-26
HmAb_Pierce2016{14} = [441 442 494 497 502:509 511 537 539 552 554 559 564 603 607 608 611 614 617:619 621 640 643 644];

%HC33-1
HmAb_Pierce2016{15} = [413 418 420];
%HC33-4
HmAb_Pierce2016{16} = [408 413 420];
% %CD81_bs
% HmAb_Pierce2016{17} = [420 421 424 427 430 436:438 440:443 523 527 529 530 535 540 613 614 616:618];



%HCV1
HmAbs_Gopal2017{1} = [413 415 418 420 422 624];
%AR1A
HmAbs_Gopal2017{2} = [417 485 494 497 502 504 506:511 514 516 519 537:540 542 544 545 547 549:552 554 559 564 565 602 603 607:611 614 617:619 621 623 624 632 638 640 642 643 644]; %559 is RB=21%
%AR1B
HmAbs_Gopal2017{3} = [494 497 504 506:509 511 537 539 544 545 547:552 554 564 607:608 610 611 614 617:619 621 640 644]; %610 is RB=21%
%AR2A
HmAbs_Gopal2017{4} = [552 607:611 614 617:619 621 623:625 628 638 640 643 644];
%AR3A
HmAbs_Gopal2017{5} = [424 425 428 429 436:438 440:442 485 494 496 497 502:509 511 516 523 525 529 530 535 537 539 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 640 642:644];
%AR3B
HmAbs_Gopal2017{6} = [424 425 427:429 432 436 437 440:442 494 496 497 502:509 511 516 518 520 523 529 530 535 537 539 540 552 554 555 556 558 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640:644];
%AR3C
HmAbs_Gopal2017{7} = [424 425 428:429 437 438 441:443 494 496 497 502:509 511 516 525 529 530 535 537 539 540 552 554 559 564 602 603 607:611 614 616:619 621 623 624 625 630 638 640 641 643 644];
%AR3D
HmAbs_Gopal2017{8} = [424 425 427:430 436:438 440:442 459 494 496 497 499 502:509 511 516 518 520 523 529 530 535 537 539 540 550:552 554 555 556 558 559 564 565 600 602 603 607:612 614 616:619 621 623 624 625 638 640:644];
 
 

%CBH-5
HmAbs_Keck2019{1} = [424:429 436:438 441:443 535 638];
%212.1.1
HmAbs_Keck2019{2} = [425:429 433 434 529 530 535];
%212.1.10
HmAbs_Keck2019{3} = [426 428 429 442];
%212.15
HmAbs_Keck2019{4} = [544 547 549 637 638 639];
%212.25
HmAbs_Keck2019{5} = [549 636 638 639];

 
% Bailey2017

%HEPC3
HmAbs_Bailey2017{1} = [425 427 428 437 499 520 530 535];
%HEPC43
HmAbs_Bailey2017{2} = [425 427 428 432 436 437 438 442 443 499 517 520 527 529 530 535 616];
%HEPC74
HmAbs_Bailey2017{3} = [425 428 436 437 530 535];
%HEPC46
HmAbs_Bailey2017{4} = [541:546 548 549 594 598 633];
%HEPC50
HmAbs_Bailey2017{5} = [543 544 545 549 594 597 598];
%HEPC98
HmAbs_Bailey2017{6} = [402 405 408];

HmAbs_all=[HmAbs_all  HmAbs_Gopal2017 HmAbs_Keck2019 HmAbs_Bailey2017 HmAb_Pierce2016];


for kk = 1:length(HmAbs_all)
    data_mean_HmAb_new{kk} = all_E(HmAbs_all{kk}-383); 
    min_data_mean_HmAb_1a(kk) =min(data_mean_HmAb_new{kk});
end

load('mean_escape_time_1b.mat')
all_E=all_mean_escape_time_1b;

for kk = 1:length(HmAbs_all)
    data_mean_HmAb_new{kk} = all_E(HmAbs_all{kk}-383); 
    min_data_mean_HmAb_1b(kk) =min(data_mean_HmAb_new{kk});
end


data = [min_data_mean_HmAb_1a  min_data_mean_HmAb_1b ];
G = [zeros(size(min_data_mean_HmAb_1a )) ones(size(min_data_mean_HmAb_1b))];



set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
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
outlier_marker_edgeColor = 'w';label_xaxis_data = {'1a',sprintf('1b')};
text_ylabel = {'Minumum escape time of',' binding resiudes of HmAbs'};
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 0.2;
savefig = 0;
savefig_name = 'escape_time_compare';
fig_width_cm = 4;
fig_height_cm = 5;
FIG=figure;
set(gcf,'renderer','Painters')

x1=0.8+0.4*(rand(length(min_data_mean_HmAb_1a),1));
x2=1.8+0.4*(rand(length(min_data_mean_HmAb_1b),1));
size_marker = 10;
% green = [27 129 62]/255;
% lightgreen = (1-green)*0.6+green;
f1=scatter(x1,min_data_mean_HmAb_1a ,'o','MarkerEdgeColor','w','MarkerFaceColor',purple,'SizeData',size_marker,'LineWidth',0.01);f1.MarkerFaceAlpha = 0.6;hold on 
f2=scatter(x2,min_data_mean_HmAb_1b,'o','MarkerEdgeColor','w','MarkerFaceColor',orange,'SizeData',size_marker,'LineWidth',0.01);f2.MarkerFaceAlpha = f1.MarkerFaceAlpha;hold on





figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);

% V=violinplot(data, [repmat("1a",1,length(min_data_mean_HmAb_1a )) repmat("1b",1,length(min_data_mean_HmAb_1b))],'EdgeColor' ,[0 0 0],'BoxColor' ,[0 0 0]);
% SizeData=10;
% ylabel({'Minimum escape time of',' binding residues of HmAbs'});
% V(1, 1).ViolinColor = purple;
% V(1, 2).ViolinColor = orange;
% V(1, 1).EdgeColor = 'None';
% 
% V(1, 1).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).MedianPlot.SizeData  =SizeData+30;
% V(1, 2).EdgeColor = 'None';
% V(1, 1).ScatterPlot.SizeData  =SizeData;
% V(1, 2).ScatterPlot.SizeData  =SizeData;

set(gca,'YTick',0:200:400)
yt = get(gca, 'YTick');
axis([xlim    0  480])
xt = get(gca, 'XTick');
hold on
plot(xt([1 2]), [1 1]*440, '-k','LineWidth',0.3)
% plot(xt([1 1]), [0.95 1]*440, '-k','LineWidth',0.3)
% plot(xt([2 2]), [0.95 1]*440, '-k','LineWidth',0.3)
P = ranksum(min_data_mean_HmAb_1b,min_data_mean_HmAb_1a,'tail','left')

ind = floor(log10(P));
P = roundn(P/10^ind,-1);
t = ['$$ P = ',num2str(P),' \times 10^{',num2str(ind),'} $$'];

% text(0.9,480,t,'interpreter','latex','FontSize',8)


FIG.Name = 'escape_time_compare';


FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 8 8.58]);
set(gcf,'Position',[6.53 6.53 4 6.27]);
% set(gcf,'Position',[6.53 6.53 5.22 3.92]);
set(gca,'TickLength',[0.035, 0.01])
set(gca,'Position',[.37 .1 .6 .88]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(gca,'TickLength',[0.02, 0.01])
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.5, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpdf');     end
