% the captions on the figure are manually adjusted, and the the circle (1a
% E2 and 1b E2) was colored using Photoshop

[header,msa] = fastaread('AL_1b_polyprotein.fasta');

 msa = cell2mat(msa');
 phi_opt = 0.88;
 all_num = [];
 all_len =[];
 
 
 all_entropy = [];
%% Core 1-191

seq = msa(:,1:191);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);
entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];








 all_entropy = [all_entropy;entropy];

%% E1 192-383

seq = msa(:,192:383);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);


entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];



 all_entropy = [all_entropy;entropy];
 
 
%% E2 384:746

seq = msa(:,384:746);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);

entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);

k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];



 all_entropy = [all_entropy;entropy];
 
 %% p7 747:809

seq = msa(:,747:809);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);

entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);

k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];



 all_entropy = [all_entropy;entropy];
 
%% NS2 810:1026

seq = msa(:,810:1026);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);

entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];




 all_entropy = [all_entropy;entropy];
 
%% NS3 1027:1657

seq = msa(:,1027:1657);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);

entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];




 all_entropy = [all_entropy;entropy];
 
%% NS4A 1658:1711

seq = msa(:,1658:1711);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);

entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];




 all_entropy = [all_entropy;entropy];
 
 %% NS4B  1712:1972

seq = msa(:,1712:1972);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);


entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];



 all_entropy = [all_entropy;entropy];
 
  %% NS5A  1973:2420

seq = msa(:,1973:2420);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);


entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];



 all_entropy = [all_entropy;entropy];
 
   %% NS5b  2420:3011

seq = msa(:,2421:end);
 all_len =[all_len;size(seq,2)];
[Profile, Symbols] = seqprofile(seq);


entropy = Profile.*log(Profile); entropy = sum(entropy(~isnan(entropy)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para = (sum(k)^2+sum(k))/2;
all_num = [all_num;num_of_para];



 all_entropy = [all_entropy;entropy];
 [header,msa] = fastaread('AL_1a_polyprotein.fasta');
%% 1a E2
msa = cell2mat(msa');

seq = msa(:,384:746);
[Profile, Symbols] = seqprofile(seq);
entropy_1a = Profile.*log(Profile); entropy_1a = sum(entropy_1a(~isnan(entropy_1a)))/size(seq,2);
k = sum(sum(Profile>0)-1); num_of_para_1a = (sum(k)^2+sum(k))/2;