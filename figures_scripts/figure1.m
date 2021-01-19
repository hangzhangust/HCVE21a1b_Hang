%% generate the correlation of the max-ent model
clear;
load('Model_1b.mat')
run startup.m
%1
[header,align] = fastaread('23890816_E2.fasta');
seq=align(2:3);
seq=cell2mat(seq');
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
end

E = (E-mean(E))/std(E);
E1_norm = E;
I = [79139.9;31818.1];
I1_norm=I;
I1_norm = (I1_norm-mean(I1_norm))./std(I1_norm);
I=[89510;38746.8];
I = log(I./I(1));
I = (I - mean(I))./(std(I));
I2_norm =I;
[r1,p1] = corr(E1_norm,I1_norm,'type','spearman');

E1_norm=E1_norm;
I1_norm=I1_norm;
E2_norm = E;
[r2,p2] = corr(E2_norm,I2_norm,'type','spearman');
%%
%2
[header,align] = fastaread('26699643_E2.fasta');
seq=align(8:11);
seq=cell2mat(seq');
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
end

E = (E-mean(E))/std(E);
E3_norm = E;
I=[282.475;1022.83;696.847;762.699];
I3_norm=I;
I3_norm = (I3_norm-mean(I3_norm))./std(I3_norm);

[r3,p3] = corr(E3_norm,I3_norm,'type','spearman');

%%
%3_1
[~,align] = fastaread('28558001.fasta');
seq=align(2:3);
seq=cell2mat(seq');
seq = seq(:,23:end);
seq = seq(:,193:end);
b_2 = seq(1,:);
b_2_1=b_2;
b_2_1(5)='V';
b_2_2=b_2;
b_2_2(12)='A';
b_2_3=b_2;
b_2_3(5)='V';
b_2_3(12)='A';
seq = [b_2;b_2_1;b_2_2;b_2_3];
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
end

E = (E-mean(E))/std(E);
E4_norm = E;
I = [1;1.11;1.41;0.828];
I4_norm=I;
I4_norm = (I4_norm-mean(I4_norm))./std(I4_norm);


[r4,p4] = corr(E4_norm,I4_norm,'type','spearman');
%%
%3_2
[header,align] = fastaread('28558001.fasta');
seq=align(2:3);
seq=cell2mat(seq');
seq = seq(:,23:end);
seq = seq(:,193:end);
b_3 = seq(2,:);
b_3_1 = seq(1,:);
b_3_1(5)='V';
b_3_2 = b_2;
b_3_2(12)='A';
b_3_3 = b_2;
b_3_1(5)='V';
b_3_3(5)='V';
b_3_3(12)='A';
seq = [b_3;b_3_1;b_3_2;b_3_3];
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i)= seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
end

E = (E-mean(E))/std(E);
E5_norm = E;
I = [0.537;0.527;0.758;0.591];
I5_norm=I;
I5_norm = (I5_norm-mean(I5_norm))./std(I5_norm);


[r5,p5] = corr(E5_norm,I5_norm,'type','spearman');

% %%
% % load('23706810_E2');
% % seq = seqE2([4 5],:);
% % seq(1,14)='W';
% % [header,align] = fastaread('new_1b_E2.fasta');
% % seq=align;
% % seq=cell2mat(seq');
% seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
% E = zeros(size(seq,1),1);
% for i =1:size(seq,1)
% % E(i)= seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
% E(i)= seq_bin(i,:)*H';
% end
% 
% E = (E-mean(E))/std(E);
% E6_norm = E;
% % I = [1.27532;2.42583];
% I = [33.07555855;1418.896504];
% 
% I6_norm=I;
% I6_norm = (I6_norm-mean(I6_norm))./std(I6_norm);
% 
% 
% [r6,p6] = corr(E6_norm,I6_norm,'type','spearman');
%%
% 23458406
[header,align] = fastaread('Con1_E2.fasta');
s=[];
seq=align;
s= [s;seq];
seq(417-383) = 'S';
s = [s;seq];
seq(417-383) = 'T';
s = [s;seq];
seq(417-383) = 'G';
s = [s;seq];
seq=align;
seq(415-383) = 'D';
s = [s;seq];
seq=s;
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i)= seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
% E(i)= seq_bin(i,:)*H';
end

E = (E-mean(E))/std(E);
E6_norm = E;
I = [157995820.3;10316776.4;3744237.867;115665920.3;7552727.877];

I6_norm=I;
I6_norm = (I6_norm-mean(I6_norm))./std(I6_norm);


[r6,p6] = corr(E6_norm,I6_norm,'type','spearman');
%% plot figures

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)

markersize = 6;
line_width = 0.25;
FIG=figure;
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
% xlabel('$${\rm E}_{\rm mutant}-{\rm E}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$${\rm f}_{\rm mutant}/{\rm f}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$$\frac{f_{mutant}}{f_{Mahoney}}$$','interpreter','latex')
hold on 
plot(E1_norm-0.2,I1_norm-0.2,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',red,'LineWidth',line_width)
plot(E2_norm+0.3,I2_norm+0.3,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',orange,'LineWidth',line_width)
plot(E3_norm,I3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',green,'LineWidth',line_width)
plot(E4_norm,I4_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',blue,'Linewidth',line_width)
plot(E5_norm,I5_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
plot(E6_norm,I6_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',pink,'LineWidth',line_width)
% plot(E7_norm,I7_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',[0 0 0],'LineWidth',line_width)
% h = legend('weigt2008','weigt2008','dunn2007','morcos2011','morcos2011');
% h.Position = [0.70 0.75 0.12 0.0488];
total_length = length(I1_norm)+length(I2_norm)+length(I3_norm)+length(I4_norm)+length(I5_norm)+length(I6_norm);

w1 = length(I1_norm)/total_length;
w2 = length(I2_norm)/total_length;
w3 = length(I3_norm)/total_length;
w4 = length(I4_norm)/total_length;
w5 = length(I5_norm)/total_length;
w6 = length(I6_norm)/total_length;
% w7 = length(I7_norm)/total_length;
W = [w1*ones(length(I1_norm),1);w2*ones(length(I2_norm),1);w3*ones(length(I3_norm),1);w4*ones(length(I4_norm),1);w5*ones(length(I5_norm),1);w6*ones(length(I6_norm),1)];

rho_weighted_average = r1*w1+r2*w2+r3*w3+r4*w4+r5*w5+r6*w6

% text(2.1,-1,sprintf('$$\\bar{r} = %.2f$$',rho_weighted_average),'interpreter','latex','FontSize',12)

P = lscov([[E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm] ones(total_length,1)],[I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],W)
% P = polyfit([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]',1);
% P = polyfit([E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm],...
%     [I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],1);
x = -2:.5:3; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1)
% h = legend('(34)','(34) M299V','(35)','(36)','(36) T292A','(37)','Location','eastoutside','Orientation','horizontal','NumColumns',1);
h = legend('(33)','(33) M299V','(34)','(35)','(35) T292A','(36)','Location','eastoutside');
% h = legend('(25)','(25) M299V','(26)','(27)','(27) T292A','(28)','Location','eastoutside','Orientation','vertical');
legend boxoff
set(h,...
    'Position',[0.680104166666667 0.331023372843418 0.30646188494189 0.467344682712138],...
    'Orientation','horizontal');
ylim([-2 3]);
xlim([-2 5])
xticks([-4:1:5])
% legend boxoff
FIG.Name = '579_model'; 
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 10 6]);
set(gca,'Position',[.11 .18 .6 .77]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.16, 0.5, 0]);
set(findobj('FontSize',8),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
try 
    print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 

catch
   
end

%% generate the correlation of the con-only model
clear;
%% conservation-only
load('Model_1b.mat')
run startup.m
Single_Mutation_Oberved = sum(msa_bin.*weight_seq,1)/sum(weight_seq);
H = Single_Mutation_Oberved;
H = -log((H)./(1-H));
%1
[~,align] = fastaread('23890816_E2.fasta');
seq=align(2:3);
seq=cell2mat(seq');

seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*H';
end

E = (E-mean(E))/std(E);
E1_norm = E;
I = [79139.9;31818.1];
I1_norm=I;
I1_norm = (I1_norm-mean(I1_norm))./std(I1_norm);
I=[89510;38746.8];
I = I;
I = (I - mean(I))./(std(I));
I2_norm =I;
[r1,p1] = corr(E1_norm,I1_norm,'type','spearman');

E1_norm=E1_norm;
I1_norm=I1_norm;
E2_norm = E;
[r2,p2] = corr(E2_norm,I2_norm,'type','spearman');
%%
%2
[header,align] = fastaread('26699643_E2.fasta');
seq=align(8:11);
seq=cell2mat(seq');
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*H';
end

E = (E-mean(E))/std(E);
E3_norm = E;
I=[282.475;1022.83;696.847;762.699];
I3_norm=I;
I3_norm = (I3_norm-mean(I3_norm))./std(I3_norm);

[r3,p3] = corr(E3_norm,I3_norm,'type','spearman');

%%
%3_1
[header,align] = fastaread('28558001.fasta');
seq=align(2:3);
seq=cell2mat(seq');
seq = seq(:,23:end);
seq = seq(:,193:end);
b_2 = seq(1,:);
b_2_1=b_2;
b_2_1(5)='V';
b_2_2=b_2;
b_2_2(12)='A';
b_2_3=b_2;
b_2_3(5)='V';
b_2_3(12)='A';
seq = [b_2;b_2_1;b_2_2;b_2_3];
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*H';
end

E = (E-mean(E))/std(E);
E4_norm = E;
I = [1;1.11;1.41;0.828];
I4_norm=I;
I4_norm = (I4_norm-mean(I4_norm))./std(I4_norm);
[r4,p4] = corr(E4_norm,I4_norm,'type','spearman');


%%
%3_2
[header,align] = fastaread('28558001.fasta');
seq=align(2:3);
seq=cell2mat(seq');
seq = seq(:,23:end);
seq = seq(:,193:end);
b_3 = seq(2,:);
b_3_1 = seq(1,:);
b_3_1(5)='V';
b_3_2 = b_2;
b_3_2(12)='A';
b_3_3 = b_2;
b_3_1(5)='V';
b_3_3(5)='V';
b_3_3(12)='A';
seq = [b_3;b_3_1;b_3_2;b_3_3];
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
E(i) = seq_bin(i,:)*H';
end

E = (E-mean(E))/std(E);
E5_norm = E;
I = [0.537;0.527;0.758;0.591];
I5_norm=I;
I5_norm = (I5_norm-mean(I5_norm))./std(I5_norm);


[r5,p5] = corr(E5_norm,I5_norm,'type','spearman');

%%
% 23458406
[header,align] = fastaread('Con1_E2.fasta');
s=[];
seq=align;
s= [s;seq];
seq(417-383) = 'S';
s = [s;seq];
seq(417-383) = 'T';
s = [s;seq];
seq(417-383) = 'G';
s = [s;seq];
seq=align;
seq(415-383) = 'D';
s = [s;seq];
seq=s;
seq_bin=Binary_Seq(seq,amino_single_combine_array,conserved);
E = zeros(size(seq,1),1);
for i =1:size(seq,1)
% E(i)= seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
E(i)= seq_bin(i,:)*H';
end

E = (E-mean(E))/std(E);
E6_norm = E;
I = [157995820.3;10316776.4;3744237.867;115665920.3;7552727.877];

I6_norm=I;
I6_norm = (I6_norm-mean(I6_norm))./std(I6_norm);


[r6,p6] = corr(E6_norm,I6_norm,'type','spearman');



%%
set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',8)
set(0,'DefaultTextFontSize',8)


line_width = 0.25;
FIG=figure;
% xlabel('Energy (normalized)')
% ylabel({'Experimental fitness';'(normalized)'})
% xlabel('$${\rm E}_{\rm mutant}-{\rm E}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$${\rm f}_{\rm mutant}/{\rm f}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$$\frac{f_{mutant}}{f_{Mahoney}}$$','interpreter','latex')
hold on 

markersize = 4;



plot(E1_norm-0.2,I1_norm-0.2,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',red,'LineWidth',line_width)
plot(E2_norm,I2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',orange,'LineWidth',line_width)
plot(E3_norm,I3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',green,'LineWidth',line_width)
plot(E4_norm,I4_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',blue,'Linewidth',line_width)
plot(E5_norm,I5_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
plot(E6_norm,I6_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',pink,'LineWidth',line_width)
% h = legend('weigt2008','weigt2008','dunn2007','morcos2011','morcos2011');
% h.Position = [0.70 0.75 0.12 0.0488];
total_length = length(I1_norm)+length(I2_norm)+length(I3_norm)+length(I4_norm)+length(I5_norm)+length(I6_norm);

w1 = length(I1_norm)/total_length;
w2 = length(I2_norm)/total_length;
w3 = length(I3_norm)/total_length;
w4 = length(I4_norm)/total_length;
w5 = length(I5_norm)/total_length;
w6 = length(I6_norm)/total_length;
rho_weighted_average = r1*w1+r2*w2+r3*w3+r4*w4+r5*w5+r6*w6
W = [w1*ones(length(I1_norm),1);w2*ones(length(I2_norm),1);w3*ones(length(I3_norm),1);w4*ones(length(I4_norm),1);w5*ones(length(I5_norm),1);w6*ones(length(I6_norm),1)];
% text(-1.9,-2.4,sprintf('$$\\bar{r} = %.2f$$',rho_weighted_average),'interpreter','latex','FontSize',8)

% P = polyfit([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]',1);
% P = polyfit([E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; I6_norm],...
%     [I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],1);
P = lscov([[E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm; E6_norm] ones(total_length,1)],[I1_norm;I2_norm;I3_norm;I4_norm;I5_norm; I6_norm],W)
% corr([E1_norm; E2_norm ;E3_norm; E4_norm; E5_norm],[I1_norm;I2_norm;I3_norm;I4_norm;I5_norm],'type','spearman');


x = -1.5:.5:2; %xaxis
y = polyval(P,x);
plot(x,y,'k--','LineWidth',1)



% h = legend('(9)','(9)','(10)','(11)','(11)');

% h = legend('(1)','(1) M491V','(2)','(3)','(3) T484A','(4)','Location','northeast');
% ylim([-2 2]);
% xlim([-2 3])
% xticks([-2:1:3])
ylim([-2.8 2]);
xlim([-2.5 2.5])
xticks([-4:1:4])
% ylim([-2 2]);
% legend boxoff

% FIG.Units = 'centimeters';
% FIG.Name = '579_con';
% set(gcf,'Position',[10 10 3.8 3.8]);
% set(gca,'Position',[.29 .2 .68 .75]);  %调整 XLABLE和YLABLE不会被切掉
% figure_FontSize=8;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.3, 0.5, 0]);
% set(get(gca,'XLabel'), 'Units', 'Normalized', 'Position', [0.47, -0.13, 0]);
% set(findobj('FontSize',8),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
% print(['C:\Users\27909\Desktop\' FIG.Name],'-depsc','-r600');


% h = legend('(9)','(9) M491V','(10)','(11)','(11) T484A');
% legend boxoff




title({'Conservation-only model'})
% ylabel({'Experimental fitness (normalized)'})
% ylim([-2 2]);
FIG.Name = '579_con';
FIG.Units = 'centimeters';
set(gcf,'Position',[10 10 3.5 2.76]);
set(gca,'Position',[.17 .17 .8 .65]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(gca,'TickDir','out')
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.14, 0.5, 0]);
set(get(gca,'title'), 'Units', 'Normalized', 'Position', [0.42, 1.1, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
set(gca,'TickLength',[0.02, 0.03])
try 
    print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 

catch

end


%% Contact Prediction (Frobenius norm)

num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);


num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
num = size(msa_aa,2);
F = zeros(num,num);
% J_MPF_BML = 2*J_MPF_BML-diag(diag(J_MPF_BML));
for i = 1:num
    Start_site_i = num_mutants_combine_array_acc_all(i)+1;
    End_site_i   = num_mutants_combine_array_acc_all(i+1);
   
    for j =1:num
        Start_site_j = num_mutants_combine_array_acc_all(j)+1;
        End_site_j   = num_mutants_combine_array_acc_all(j+1); 
        F(i,j)=F(i,j)+sum(sum((J_MPF_BML(Start_site_i:End_site_i,Start_site_j:End_site_j).^2)));
    
    end
end
% F = (F+F')/2;
F = F.^0.5;
F_new=F;
avg=0;
for i = 1:num  
    for j =1:num
        if i ~=j
            avg = avg+F(i,j);
        end
    end
end

% avg = mean(mean(F));
avg = avg/(num^2-num);



for i =1:num
    for j =1:num
        
          avg_i = mean(F(:,j))-F(j,j)/num;
          avg_j = mean(F(i,:))-F(i,i)/num;
%         avg_i = mean(F(:,j));
%         avg_j = mean(F(i,:));
        F_new(i,j) = F(i,j) - avg_i*avg_j/avg;
    end
end
F = F_new;
F =triu(F,1);


Nmax = num*(num-1)/2; % get Nmax biggest entries
[ ~, Ind ] = sort(F(:),1,'descend');

[ ind_row, ind_col ] = ind2sub(size(F),Ind(1:Nmax)); % fetch indices
load('6mei.mat')
indx=indx-383;
indy=indy-383;

% load('I-TASSER\model5.mat')
ind_non_conserve = setdiff(1:363,conserved);
ind_col = ind_non_conserve(ind_col);
ind_row = ind_non_conserve(ind_row);

ia1 = ind_row<=(647-383) & ind_row>= (410-383);
ia2 = ind_col<=(647-383) & ind_col>= (410-383);
ia = ia1 & ia2;
ind_col =ind_col(ia);
ind_row =ind_row(ia);


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




FIG=figure;
plot(TPR,'Color','k','LineWidth',1,'Color',[0.6 0.6 0.6])
hold on


set(gca, 'XScale', 'log')


xlabel('Top x pairs')
ylabel({'Precision'})
FIG.Units = 'centimeters';
FIG.Name = 'Ours';

box off
xlim([1 1e3])

xticks([1 10 100 1000])
xticklabels({'10^{0}','10^{1}','10^{2}','10^{3}'})
FIG.Units = 'centimeters';
% set(gcf,'Position',[10 10 7.84 6]);
% set(gca,'Position',[.15 .2 .8 .74]);  %调整 XLABLE和YLABLE不会被切掉
set(gcf,'Position',[10 10 5.3 6]);
set(gca,'Position',[.24 .18 .7 .77]);  %调整 XLABLE和YLABLE不会被切掉
figure_FontSize=8;
set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
set(gca,'TickDir','out')
set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
set(get(gca,'YLabel'), 'Units', 'Normalized', 'Position', [-0.26, 0.5, 0]);
set(findobj('FontSize',10),'FontSize',figure_FontSize);
set(gca,'TickLength',[0.02, 0.03])
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
try 
    print(['C:\Users\27909\Desktop\' FIG.Name],'-dpng','-r600'); 

catch

end