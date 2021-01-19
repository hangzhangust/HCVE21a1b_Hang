%% all 35 antibodies 
%% 1a
run startup.m
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

HmAbs_all=[HmAbs_all  HmAbs_Bailey2017 HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016];


for kk = 1:length(HmAbs_all)
    data_mean_HmAb_new{kk} = all_E(HmAbs_all{kk}-383); 
    min_data_mean_HmAb_new(kk) =min(data_mean_HmAb_new{kk});
end

ylim_min = 0;
ylim_max = 400;




FIG = figure;
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
% common 
% 
% xbars = [0.5 1.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% 
% xbars = [11.5 13.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% xbars = [17.5 18.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% 
% xbars = [25.5 26.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [28.5 30.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [33.5 34.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)


%1a
% xbars = [0.5 1.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [1 1 1], ...
%     'EdgeColor','w','LineWidth',0.1)
% box off

%1a
% xbars = [0.5 1.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% 
% % xbars = [10.5 13.5];
% 
% % xbars = [11.5 14.5];
% % patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
% %     'EdgeColor','w','LineWidth',0.1)
% xbars = [6.5 7.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [2.5 5.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [17.5 19.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [25.5 26.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [28.5 31.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [33.5 34.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [0 35.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[1 1 1], ...
%     'EdgeColor','w','LineWidth',0.1)








% %1b
% xbars = [0.5 1.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% 
% xbars = [10.5 13.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% xbars = [17.5 18.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% 
% xbars = [25.5 26.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [28.5 30.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [33.5 34.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% blue = color_scheme_npg(2,:);
% 
% red = color_scheme_npg(8,:);

h1 = barh(1:35,[zeros(1,1) min_data_mean_HmAb_new(2) zeros(1,3) min_data_mean_HmAb_new(6) zeros(1,29)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
hold on
barh(1:35,[min_data_mean_HmAb_new(1) zeros(1,1) min_data_mean_HmAb_new(3:5) zeros(1,1)  zeros(1,29)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);

h2=barh(1:35, [zeros(1,6) min_data_mean_HmAb_new(7) zeros(1,7) zeros(1,21)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);
barh(1:35, [zeros(1,6) zeros(1,1) min_data_mean_HmAb_new(8:14) zeros(1,21)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);

h3 = barh(1:35, [zeros(1,14) zeros(1,3) min_data_mean_HmAb_new(18:19) zeros(1,16)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);
barh(1:35, [zeros(1,14) min_data_mean_HmAb_new(15:17) zeros(1,2) zeros(1,16)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
% h2 = barh(1:35, [zeros(1,6) min_data_mean_HmAb_new(7:14) zeros(1,21)], 0.4, 'FaceColor',color_scheme_npg(7,:),'FaceAlpha',0.7);
% h3 = barh(1:35, [zeros(1,14) min_data_mean_HmAb_new(15:19) zeros(1,16)], 0.4, 'FaceColor',color_scheme_npg(8,:),'FaceAlpha',0.7);
h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:25)  zeros(1,1) min_data_mean_HmAb_new(27:28) zeros(1,3) min_data_mean_HmAb_new(32:33) zeros(1,1)  min_data_mean_HmAb_new(35)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
barh(1:35, [zeros(1,19) zeros(1,6)  min_data_mean_HmAb_new(26) zeros(1,2) min_data_mean_HmAb_new(29:31) zeros(1,2) min_data_mean_HmAb_new(34)  zeros(1,1)], 0.4, 'FaceColor',red,'FaceAlpha',0.7); 
% h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:35) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);



% h1 = barh(1:35,[min_data_mean_HmAb_new(1:8) zeros(1,27)], 0.4, 'FaceColor',color_scheme_npg(6,:),'FaceAlpha',0.7);
% hold on
% h2 = barh(1:35, [zeros(1,8) min_data_mean_HmAb_new(9:13) zeros(1,22)], 0.4, 'FaceColor',color_scheme_npg(7,:),'FaceAlpha',0.7);
% h3 = barh(1:35, [zeros(1,13) min_data_mean_HmAb_new(14:19) zeros(1,16)], 0.4, 'FaceColor',color_scheme_npg(8,:),'FaceAlpha',0.7);
% h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:35) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
ylim([ylim_min ylim_max])
set(gca,'TickDir','out')

ylim([0 35.5])



% yt = get(gca, 'YTick');
% xt = get(gca, 'XTick');
% plot([70 70], [0 35.5], '--k','LineWidth',0.5)
% text(150,36,'Subtype 1b','FontSize',8)



%1a

box off
yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
plot([100 100], [0 35.5], '--k','LineWidth',0.5)
% text(120,36,'Subtype 1a','FontSize',8)
% xtick =set(gca,'YTick',1:35,'YTickLabel',...
%     {'\bf HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
%     'AR3C','AR3D','CBH-5','212.1.1',' 212.10','\bf 212.15','\bf 212.25','HEPC3',...
%     'HEPC43',' HEPC74',' HEPC46','\bf HEPC50','HEPC98','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
%     '\bf HC-1','HC-11',' A27','\bf CBH-23','\bf CBH-7',' HC84-20','HC84-24','HC84-26',...
%     '\bf HC33-1','HC33-4'},'FontName','Arial');

ytickangle(315)

xticks([0 100  400])
xtick =set(gca,'YTick',[],'YTickLabel',...
    {},'FontName','Arial');
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');

% h = legend([h1(1), h2(1) ,h3(1) ,h4(1)],'Gopal et al.','Keck et al.','Bailey et al.','Pierce et al.' ,'Location','southoutside','Orientation','horizontal');

%1b
% xlabel('Minimum escape time, t^{min}_{e}')


% set(gca,'TickLabelInterpreter','latex');
% xtick =set(gca,'YTick',[],'YTickLabel',...
%     {},'FontName','Arial');


text(0.3,1.02,'1a','Units','normalized')



FIG.Units = 'centimeters';
FIG.Name = 'RBall_1a';
set(gcf,'Position',[5 5 3 15.42]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.34 .035 .58 .93]);  %调整 XLABLE和YLABLE不会被切掉

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
% 
% label_xaxis_data = {'HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
%     'AR3C','AR3D','CBH-5','212.1.1','212.10','212.15','212.25','HEPC3',...
%     'HEPC43','HEPC74','HEPC46','HEPC50','HEPC98','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
%     'HC-1','HC-11','A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
%     'HC33-1','HC33-4'};
% 
% 
% 
% label_yaxis_data = {80,100,120,140};
% 
% 
% 
% map = zeros(length(label_yaxis_data),length(label_xaxis_data));
% for i =1:length(label_yaxis_data)
% for kk = 1:length(HmAbs_all)
%     data = all_E(HmAbs_all{kk}-383);
% %     map(i,kk) = length(find(data<= prctile(all_E,i*5)));
%     map(i,kk) = length(find(data<=label_yaxis_data{i}));
% end
% end
% 
% % map  = flip(map,1);
% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',8)
% set(0,'DefaultTextFontSize',8)
% 
% FIG=figure;
% % subplot(2,1,2)
% [grad,im]=colorGradient(	hex2rgb('#FCFAF2'),hex2rgb('#4A225D'),128);
% heatmap(label_xaxis_data,label_yaxis_data,map,'GridVisible','off','Colormap',grad,'ColorLimits',[0 max(max(map))],'ColorbarVisible','off','FontName','Arial','FontSize',8);
% % text('Units', 'Normalized', 'Position', [-0.09,0.5],'String','Escape time threshold, \tau','Rotation',90)
% % ,'FontName','Arial','FontSize',8,'FontWeight','Bold'
% 
% FIG.Units = 'centimeters';
% FIG.Name = 'RBall_thre_1b';
% set(gcf,'Position',[10 10 18 4]);
% % ylabel('Escape time threshold, \tau')
% set(gca,'Position',[.05 .3 .9 .6]);  %调整 XLABLE和YLABLE不会被切掉
% 
% 
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');



%% 1b
load('mean_escape_time_1b.mat')
all_E=all_mean_escape_time_1b;
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

HmAbs_all=[HmAbs_all HmAbs_Bailey2017  HmAbs_Gopal2017 HmAbs_Keck2019 HmAb_Pierce2016];


for kk = 1:length(HmAbs_all)
    data_mean_HmAb_new{kk} = all_E(HmAbs_all{kk}-383); 
    min_data_mean_HmAb_new(kk) =min(data_mean_HmAb_new{kk});
end

ylim_min = 0;
ylim_max = 400;




FIG = figure;
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)
% common 

% xbars = [0.5 1.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% 
% xbars = [11.5 13.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% xbars = [17.5 18.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% 
% xbars = [25.5 26.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [28.5 30.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [33.5 34.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 


%1b
% xbars = [4.5 5.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% xbars = [6.5 7.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% % xbars = [10.5 13.5];
% % patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)], [0.85 0.85 0.85], ...
% %     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% xbars = [16.5 19.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% 
% xbars = [25.5 26.5];
% patch( [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [xbars(1) xbars(1), xbars(2) xbars(2)],[0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [28.5 30.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% xbars = [33.5 34.5];
% patch([ylim_min+.05 ylim_max ylim_max ylim_min+.05],[xbars(1) xbars(1), xbars(2) xbars(2)],  [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% blue = color_scheme_npg(2,:);

% red = color_scheme_npg(8,:);

h1 = barh(1:35,[min_data_mean_HmAb_new(1:4) zeros(1,1) min_data_mean_HmAb_new(6) zeros(1,29)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
hold on
barh(1:35,[zeros(1,4) min_data_mean_HmAb_new(5) zeros(1,1) zeros(1,29)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);
% h1 = barh(1:35,[min_data_mean_HmAb_new(1:6) zeros(1,29)], 0.4, 'FaceColor',color_scheme_npg(6,:),'FaceAlpha',0.7);
% hold on
h2 = barh(1:35, [zeros(1,6)  zeros(1,1) min_data_mean_HmAb_new(8:14) zeros(1,21)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
barh(1:35, [zeros(1,6)  min_data_mean_HmAb_new(7)  zeros(1,7) zeros(1,21)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);


h3 = barh(1:35, [zeros(1,14) min_data_mean_HmAb_new(15:16) zeros(1,3) zeros(1,16)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
barh(1:35, [zeros(1,14) zeros(1,2) min_data_mean_HmAb_new(17:19) zeros(1,16)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);


h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:25) zeros(1,1)  min_data_mean_HmAb_new(27:28) zeros(1,2) min_data_mean_HmAb_new(31:33) zeros(1,1) min_data_mean_HmAb_new(35) zeros(1,0)], 0.4, 'FaceColor',blue,'FaceAlpha',0.7);
barh(1:35, [zeros(1,19)  zeros(1,6) min_data_mean_HmAb_new(26) zeros(1,2)  min_data_mean_HmAb_new(29:30) zeros(1,3) min_data_mean_HmAb_new(34) zeros(1,1) zeros(1,0)], 0.4, 'FaceColor',red,'FaceAlpha',0.7);


% h2 = barh(1:35, [zeros(1,6) min_data_mean_HmAb_new(7:14) zeros(1,21)], 0.4, 'FaceColor',color_scheme_npg(7,:),'FaceAlpha',0.7);
% h3 = barh(1:35, [zeros(1,14) min_data_mean_HmAb_new(15:19) zeros(1,16)], 0.4, 'FaceColor',color_scheme_npg(8,:),'FaceAlpha',0.7);
% h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:35) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);



% h1 = barh(1:35,[min_data_mean_HmAb_new(1:8) zeros(1,27)], 0.4, 'FaceColor',color_scheme_npg(6,:),'FaceAlpha',0.7);
% hold on
% h2 = barh(1:35, [zeros(1,8) min_data_mean_HmAb_new(9:13) zeros(1,22)], 0.4, 'FaceColor',color_scheme_npg(7,:),'FaceAlpha',0.7);
% h3 = barh(1:35, [zeros(1,13) min_data_mean_HmAb_new(14:19) zeros(1,16)], 0.4, 'FaceColor',color_scheme_npg(8,:),'FaceAlpha',0.7);
% h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:35) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
ylim([ylim_min ylim_max])
set(gca,'TickDir','out')

ylim([0 35.5])



yt = get(gca, 'YTick');
xt = get(gca, 'XTick');
plot([80 80], [0 35.5], '--k','LineWidth',0.5)
% text(150,36,'Subtype 1b','FontSize',8)

xtick =set(gca,'YTick',1:35,'YTickLabel',...
    {'HEPC3',...
    'HEPC43',' HEPC74',' HEPC46','\bf HEPC50','HEPC98','\bf HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
    'AR3C','AR3D','CBH-5','212.1.1',' 212.10','\bf 212.15','\bf 212.25','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    '\bf HC-1','HC-11',' A27','\bf CBH-23','\bf CBH-7',' HC84-20','HC84-24','HC84-26',...
    '\bf HC33-1','HC33-4'},'FontName','Arial');


names = {'HEPC3',...
    'HEPC43',' HEPC74',' HEPC46','\bf HEPC50','HEPC98','\bf HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
    'AR3C','AR3D','CBH-5','212.1.1',' 212.10','\bf 212.15','\bf 212.25','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    '\bf HC-1','HC-11',' A27','\bf CBH-23','\bf CBH-7',' HC84-20','HC84-24','HC84-26',...
    '\bf HC33-1','HC33-4'};

% min_data_mean_HmAb_new(min_data_mean_HmAb_new>=80)

box off
ytickangle(315)

% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');

% h = legend([h1(1), h2(1) ,h3(1) ,h4(1)],'Gopal et al.','Keck et al.','Bailey et al.','Pierce et al.' ,'Location','southoutside','Orientation','horizontal');

%1b
% xlabel('Minimum escape time, t^{min}_{e}')


% set(gca,'TickLabelInterpreter','latex');
% xtick =set(gca,'YTick',[],'YTickLabel',...
%     {},'FontName','Arial');
xlim([0 400])


xticks([0 80  400])


text(0.3,1.02,'1b','Units','normalized')
FIG.Units = 'centimeters';
FIG.Name = 'RBall_1b';
% set(gcf,'Position',[5 5 3 15.42]);
% % ylabel('Escape time threshold, \tau')
% set(gca,'Position',[.34 .035 .58 .93]);  %调整 XLABLE和YLABLE不会被切掉
set(gcf,'Position',[5 5 2 15.42]);
% ylabel('Escape time threshold, \tau')
set(gca,'Position',[.34 .035 .58 .93]);  %调整 XLABLE和YLABLE不会被切掉
% set(FIG, 'DefaultTextFontSize', 8);

% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% label_xaxis_data = {'HCV1','AR1A','AR1B','AR2A','AR3A','AR3B',...
%     'AR3C','AR3D','CBH-5','212.1.1','212.10','212.15','212.25','HEPC3',...
%     'HEPC43','HEPC74','HEPC46','HEPC50','HEPC98','CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
%     'HC-1','HC-11','A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
%     'HC33-1','HC33-4'};
% 
% 
% % 
% label_yaxis_data = {50,70,100,120};
% 
% 
% 
% 
% 
% map = zeros(length(label_yaxis_data),length(label_xaxis_data));
% for i =1:length(label_yaxis_data)
% for kk = 1:length(HmAbs_all)
%     data = all_E(HmAbs_all{kk}-383);
% %     map(i,kk) = length(find(data<= prctile(all_E,i*5)));
%     map(i,kk) = length(find(data<=label_yaxis_data{i}));
% end
% end
% 
% % map  = flip(map,1);
% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',8)
% set(0,'DefaultTextFontSize',8)
% 
% FIG=figure;
% % subplot(2,1,2)
% [grad,im]=colorGradient(	hex2rgb('#FCFAF2'),hex2rgb('#4A225D'),128);
% heatmap(label_xaxis_data,label_yaxis_data,map,'GridVisible','off','Colormap',grad,'ColorLimits',[0 max(max(map))],'ColorbarVisible','off','FontName','Arial','FontSize',8);
% % text('Units', 'Normalized', 'Position', [-0.09,0.5],'String','Escape time threshold, \tau','Rotation',90)
% % ,'FontName','Arial','FontSize',8,'FontWeight','Bold'
% 
% FIG.Units = 'centimeters';
% FIG.Name = 'RBall_thre_1b';
% set(gcf,'Position',[10 10 18 4]);
% % ylabel('Escape time threshold, \tau')
% set(gca,'Position',[.05 .3 .9 .6]);  %调整 XLABLE和YLABLE不会被切掉
% 
% 
% print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600');
% % print(FIG, '-depsc2', '-r600', '-tiff', '-loose', ['C:\Users\27909\Desktop\' FIG.Name]);
% 
% %% generate the legend and x-axis
% FIG=figure;
% 
% h1 = barh(1:35,[min_data_mean_HmAb_new(1:8) zeros(1,27)], 0.4, 'FaceColor',color_scheme_npg(6,:),'FaceAlpha',0.7);
% hold on
% h2 = barh(1:35, [zeros(1,8) min_data_mean_HmAb_new(9:13) zeros(1,22)], 0.4, 'FaceColor',color_scheme_npg(7,:),'FaceAlpha',0.7);
% h3 = barh(1:35, [zeros(1,13) min_data_mean_HmAb_new(14:19) zeros(1,16)], 0.4, 'FaceColor',color_scheme_npg(8,:),'FaceAlpha',0.7);
% h4 = barh(1:35, [zeros(1,19) min_data_mean_HmAb_new(20:35) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
% 
% 
% h = legend([h4(1), h3(1) ,h2(1) ,h1(1)],'(50)','(51)','(52)','(53)' ,'Location','southoutside','Orientation','horizontal');
% ylim([0 35.5])
% xlim([0 400])
% FIG.Units = 'centimeters';
% FIG.Name = 'RBall_1b';
% set(gcf,'Position',[10 10 20 18]);
% % ylabel('Escape time threshold, \tau')
% set(gca,'Position',[.23 .1 .7 .85]);  %调整 XLABLE和YLABLE不会被切掉
% ylabel('Minimum escape time, t^{min}_{e}')
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
% % h = legend([h1(1), h2(1) ,h3(1) ,h4(1)],'Gopal et al.','Keck et al.','Bailey et al.','Pierce et al.' ,'Location','southoutside','Orientation','horizontal');
% legend boxoff
% 
% try print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); catch      print(['C:\Users\hzhangbr\Desktop\' FIG.Name],'-dpng','-r600');     end
% 
% %% boxplot of the t^min_e