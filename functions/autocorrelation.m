poolobj = gcp('nocreate');
delete(poolobj);
% cluster = parallel.cluster.Generic;
% pobj = parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
pobj = parpool(4);
% load('data_fitnessCosts_E2_99900.mat')
load('Only_confidence_best.mat')
N=1000000;
tic;
D=30;
msa_bin_unique = unique(msa_bin,'rows');
samples=randsample(size(msa_bin_unique,1),N,true);
samples = msa_bin_unique(samples,:);

% generate samples
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);


num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];

positions = size(num_mutants_combine_array,2);   
% x0=samples;
for i =1:N
    xnew =samples(i,:);
    pos = randi(positions,D,1);
    for j=1:D
    index = randi([0,num_mutants_combine_array(1,j)]);
    Start_site = num_mutants_combine_array_acc_all(j)+1;
    End_site = num_mutants_combine_array_acc_all(j+1);
        if(index==0)          
            xnew(Start_site:End_site)=0; 
              
        else
            xnew(Start_site:End_site)=0;        
            xnew(Start_site+index-1)=1;
        end
    end
    samples(i,:)=xnew;
end


steps=50;
all_Energy = zeros(N,steps+1);
%random walk
parfor i =1:N
    tmp_energy = zeros(1,steps+1);
    x = samples(i,:);
    tmp_energy(1) = x*J_MPF_BML*x';
    energy = tmp_energy(1);
    for j=1:steps
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

        else
            xnew(Start_site:End_site)=0;
            xnew(Start_site+index-1)=1;
            energy_new = energy_new+2*xnew*J_MPF_BML(:,Start_site+index-1)-J_MPF_BML(Start_site+index-1,Start_site+index-1); 
        end
     tmp_energy(j+1)=energy_new ;
     x=xnew;
     energy =energy_new;
    end
    all_Energy(i,:) = tmp_energy;
end
save(sprintf('Engergy_1b_D_%d.mat',...
        D),...
        'all_Energy');
toc;