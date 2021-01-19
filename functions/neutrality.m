poolobj = gcp('nocreate');
delete(poolobj);
% cluster = parallel.cluster.Generic;
% pobj = parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
pobj = parpool(4);
% load data
% load('data_fitnessCosts_E2_99900.mat')
load('Only_confidence_best.mat')
% generate samples
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);


num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];

positions = size(num_mutants_combine_array,2);   
N=1000000;
tic;

% load('weight_isolate_best.mat')

samples=randsample(size(msa_aa,1),N,true);

% generate samples

positions = size(num_mutants_combine_array,2);  

L=500;
epsilon = 10^-5:10^-5:10^-3;
all_d = zeros(length(epsilon),1);
ind=0;
for e = epsilon
tmp_d = zeros(1,N);
parfor i =1:N
    d = 0;
    x0 = msa_bin(samples(i),:);
    x = x0;
    energy = x*J_MPF_BML*x';
    for j=1:L
    xnew =x;
    energy_new = energy;
    pos = randi(positions);
    index = randi([0,num_mutants_combine_array(1,pos)]);
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    one = find(xnew(Start_site:End_site)>0);
    if(~isempty(one))
        energy_new = energy_new-2*xnew*J_MPF_BML(:,Start_site+one-1)+J_MPF_BML(Start_site+one-1,Start_site+one-1);
    end
        if(index==0)          
                xnew(Start_site:End_site)=0; 
            if all(xnew==x)
                j=j+1;
                continue;
            end
        else
            xnew(Start_site:End_site)=0;
            xnew(Start_site+index-1)=1;
            
            if all(xnew==x)
                j=j+1;
                continue;
            end
            
            energy_new = energy_new+2*xnew*J_MPF_BML(:,Start_site+index-1)-J_MPF_BML(Start_site+index-1,Start_site+index-1); 
        end
        if abs(energy-energy_new)<=e
        x=xnew;
        energy =energy_new;
        end
        dis=abs(sum(x0-x));
        if d<dis
            d=dis;
        end
    end
    tmp_d(i) = d;
end
ind=ind+1;
all_d(ind)=mean(tmp_d);
end
save('distance_1b_L_500.mat','all_d')
toc;