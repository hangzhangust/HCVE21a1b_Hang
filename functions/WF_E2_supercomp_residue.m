function WF_E2_supercomp_residue(imm_pres_site,strain_indx_1,strain_indx_2,ITER,Ng,b,beta_param,num_of_cores)


%
%  imm_pres_site=[415;416;417;422;424;431;433;435;438;442;444;446;447;453;456;458;461;466;475;478;482;501;524;528;531;533;538;557;558;560;580;610;636;655;713];
% after modifying imm_pres_site=[31,32,33,38,39,44,46,48,51,54,56,58,59,65,68,70,73,78,87,90,94,110,128,132,134,136,141,158,159,161,181,201,222,240,288]


% Inputs:
% imm_pres_site:    The E2 site at which immune pressure will be applied
% strain_indx_1:    The T/F strain index i in which there is WT at imm_pres_site
% strain_indx_1:    The T/F strain index k in which there is WT at imm_pres_site
%                   Simulation is run for T/F strain index i to k
% ITER:             No. of iterations to run evolutionary sim starting from
%                   the same T/F strain
% Ng:               No. of max generations to run the code
% b:                Parameter b in evol. fitness simulation that models
%                   immune pressure
% beta_param:       Parameter beta in evol. fitness eq. that models the
%                   mapping from fitness to probability of survival

%imm_pres_site = str2double(imm_pres_site)
%ITER = str2double(ITER)
%Ng = str2double(Ng)
%b = str2double(b)
%beta_param = str2double(beta_param)
%strain_indx_1 = str2double(strain_indx_1)
%strain_indx_2 = str2double(strain_indx_2)

% if maxNumCompThreads < 16
%     error('Error!No. of max. available threads on this node = %d, which is less than 16!',maxNumCompThreads)
% end
% 
% indexNumWorkers = 0;
% poolobj = gcp('nocreate');
% delete(poolobj);
% 
% poolobj = parpool(maxNumCompThreads);
% 
% while poolobj.NumWorkers < 16
%     
%     delete(poolobj);
%     
%     poolobj = parpool(maxNumCompThreads);
%     
%     poolobj.NumWorkers
%     
%     indexNumWorkers = indexNumWorkers + 1;
%     
%     if indexNumWorkers > 3
%        
%         error('Error: Assigned no. of threads on this node not more than 16!')
%         
%     end
%     
% end



poolobj = gcp('nocreate');
delete(poolobj);
% cluster = parallel.cluster.Generic;
% pobj = parpool('local', str2num(getenv('SLURM_CPUS_PER_TASK')));
pobj = parpool(num_of_cores);
% load('E1E2_1a.mat', 'J_MPF_BML')
% % load('E1E2_1a.mat', 'msa_aa')
% load('E1E2_1a.mat', 'conserved')
% load('replaced_gaps_E1E2_1a.mat', 'msa')
% load('E1E2_1a.mat', 'amino_single_combine_array')
% load('E1E2_1a.mat', 'num_mutants_combine_array')
% load('consensus_seq.mat', 'Cseq_aa')


load Model_1b.mat J_MPF_BML;
load replaced_gaps_msa_aa_1b replaced_gaps_msa_aa;
load Model_1b num_mutants_combine_array;
load Model_1b amino_single_combine_array;



% imm_pres_site:    1 to 320 [can also input multiple sites]
% strain_num:       strain number to initialize the pop (select index of random strains generated using randperm with rng(0)

H = convertRayFLParamsToJohnParamsFormat(J_MPF_BML);
msa = replaced_gaps_msa_aa;




%1a
% load('data_fitnessCosts_E2_99900.mat', 'J_mat')
% load('data_fitnessCosts_E2_99900.mat', 'amino_single_combine_array')
% load('data_fitnessCosts_E2_99900.mat', 'msa_aa')
% load('data_fitnessCosts_E2_99900.mat', 'ind_conserve')
% load('data_fitnessCosts_E2_99900.mat', 'num_mutants_combine_array')

% 
% J_mat = reshape(J_mat,[624 624]);
% H  = convertRayFLParamsToJohnParamsFormat(J_mat);
% msa = msa_aa;
% msa(:,ind_conserve) =[];
% msa(:,conserved)=[];
% 
% 
phi_curr = num_mutants_combine_array;
mutant_order = amino_single_combine_array;

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end

phi_cum = [0 phi_cumulative];

%
% msa_aa_ex = msa_bin;


% Cseq issue
Cseq_aa = seqconsensus(msa);
for kk = 1:size(msa,2)
    Cseq_mut_order(kk) = mutant_order{kk}(1);
end
% find(Cseq_mut_order~=Cseq_aa)%different because of patient weights
% %thus we use Cseq_mut_order as the consensus seq
% save Cseq_mut_order_1a Cseq_mut_order_1a
% load Cseq_mut_order Cseq_mut_order_1b
% Cseq_aa(conserved)=[];

% Modified the genetic code 23: same as standard (code = 1) except changes in code.B, code.X, and code.stop
msa_nt = aa2nt_ahmed(msa, 23);
% save msa_nt msa_nt
% load msa_nt msa_nt

func_working = ~sum(sum(msa~=nt2aa_ahmed(msa_nt, 23)));
if func_working == 1
    fprintf('\n---Modified aa2nt and nt2aa functions working as expected---\n')
end

%Parameters

N = 2000;           %Number of protein sequences in each generation
M = size(msa_nt,2); %Number of sites in the protein
mu = 1e-4;          %Mutation rate per genome per replication

%Ng = 500; %4e3;     %Number of generations the simulation will be run

imm_pres_gen = 0;
% b = 10;%12;
%beta_param = 0.1;
freq_escape = 0.5;
%ITER = 24;
%imm_pres_site = 1;
%strain_num = 1;

Cseq_nt = aa2nt_ahmed(Cseq_aa, 23);


imm_pres_site_bin = [];
for kk = 1:length(imm_pres_site)
    imm_pres_site_bin = [imm_pres_site_bin phi_cum(imm_pres_site(kk))+1:phi_cum(imm_pres_site(kk)+1)];
end
imm_pres_site_bin


%% find seqs with WT at a particular position

msa_unique = unique(msa,'rows');
msa_nt_unique = aa2nt_ahmed(msa_unique,23);

for kk = 1:size(msa_unique,2)
    seqs_wt{kk} = find(msa_unique(:,kk)==mutant_order{kk}(1)); %k-th cell contains unique MSA sequence indices with WT at k-th position
    no_seqs_wt(kk) = length(seqs_wt{kk});
end

rng(0) %initializing the random number generator

%selecting 100 random MSA sequences with WT at imm_pres_site

rand_perm = randperm(no_seqs_wt(imm_pres_site));
indx_strains = seqs_wt{imm_pres_site}(rand_perm(strain_indx_1:strain_indx_2));

strains_indices = strain_indx_1:strain_indx_2;

for strain = 1:length(indx_strains)
    
    %%
    indx_seq = indx_strains(strain);
    strain_index = strains_indices(strain);
    
    init_pop = repmat(msa_nt_unique(indx_seq,:),N,1);
    
    init_pop_aa = nt2aa_ahmed(init_pop, 23);
    % [init_pop_aa_unique,count_init_pop_aa] = unique_seqs(init_pop_aa);
    
    bin_matrix = cell(1,9);
    for aba=1:9
        temp_matrix = fliplr(eye(aba));
        temp2_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
        bin_matrix{aba} =temp2_matrix;
    end
    
    init_pop_bin = repmat(convertAAseq2Bin_ahmed(mutant_order,bin_matrix,init_pop_aa(kk,:)),N,1);
    energy_init_pop_bin = repmat(diag(init_pop_bin(kk,:)*triu(-H)*init_pop_bin(kk,:).').',N,1);
    
    % init_pop_bin_unique = zeros(size(init_pop_aa_unique,1),length(H));
    % energy_init_pop_bin_unique = zeros(size(init_pop_aa_unique,1),1);
    % energy_init_pop_bin = [];
    % for kk = 1:size(init_pop_aa_unique,1)
    %     init_pop_bin_unique(kk,:) = convertAAseq2Bin_ahmed(mutant_order,bin_matrix,init_pop_aa_unique(kk,:));
    %     energy_init_pop_bin_unique(kk) = diag(init_pop_bin_unique(kk,:)*triu(-H)*init_pop_bin_unique(kk,:).').';
    %     energy_init_pop_bin = [energy_init_pop_bin; repmat(energy_init_pop_bin_unique(kk),count_init_pop_aa,1)];
    % end
    
    parfor iter = 1:ITER
        
        %     energy_pop = cell(1,Ng+1);
        Ebar = zeros(1,Ng+1);
        Ebar(1) = mean(energy_init_pop_bin);
        
        %     Ebar_without_imm_pres = zeros(1,Ng+1);
        %     Ebar_without_imm_pres(1) = mean(energy_init_pop_bin);
        
        freq_imm_pres_site = zeros(Ng+1,1);
        freq_imm_pres_site(1,:) = 0;
        
        pop = init_pop;
        mut_pop_aa = init_pop_aa;
        mut_pop_bin = init_pop_bin;
        
        bin_matrix = cell(1,9);
        for aba=1:9
            temp_matrix = fliplr(eye(aba));
            temp2_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
            bin_matrix{aba} =temp2_matrix;
        end
        
        
        for kk = 1:Ng
            
            tic
            
            %%%%%%%%%%%%% MUTATION AT NT LEVEL %%%%%%%%%%%%%%%%%%
            mut_step = binornd(1,mu,N,M);
            
            %selecting sequences in which mutation is to be made
            seqs_with_mut = find(sum(mut_step,2).'~=0);
            
            %making nt mutations in the selected sequences randomly from ACGT
            mut_pop = pop;
            for mm = 1:length(seqs_with_mut)
                
                seq_sel = seqs_with_mut(mm);
                pos_mut = find(mut_step(seq_sel,:)~=0); %there can be multiple mutations per sequence; hence, next loop!
                
                for nn = 1:length(pos_mut)
                    
                    %nt_curr = pop(seq_sel,pos_mut(nn));
                    
                    %forming pool of nt available at this position, i.e., setdiff('ACGT',nt_curr)
                    pool_nt = setdiff('ACGT',pop(seq_sel,pos_mut(nn)));
                    %nt_new = pool_nt(randi(3));
                    
                    %selecting a random nt from the available pool
                    mut_pop(seq_sel,pos_mut(nn)) = pool_nt(randi(3));
                    
                    %finding codon indices of the mutated nucleotide
                    codon_indices = codon_indices_from_nt_indx(pos_mut(nn));
                    
                    %obtaining aa seq directly for this mutation (instead of
                    %doing it for the whole sequence)
                    mut_pop_aa(seq_sel,ceil(pos_mut(nn)/3)) = nt2aa_ahmed(mut_pop(seq_sel,codon_indices),23);
                    
                    aa_at_mut_site = mut_pop_aa(seq_sel,ceil(pos_mut(nn)/3))
                    aa_mut_site = ceil(pos_mut(nn)/3)
                    if isempty(find(mutant_order{aa_mut_site}==aa_at_mut_site, 1))
                        mut_pop_bin(seq_sel,phi_cum(aa_mut_site)+1:phi_cum(aa_mut_site+1))=...
                            bin_matrix{phi_curr(aa_mut_site)}(end,:);
                    else
                        mut_pop_bin(seq_sel,phi_cum(aa_mut_site)+1:phi_cum(aa_mut_site+1))=...
                            bin_matrix{phi_curr(aa_mut_site)}(find(mutant_order{aa_mut_site}==aa_at_mut_site, 1),:);
                    end
                    
                end
                
            end
            
            %         mut_pop_aa = nt2aa_ahmed(mut_pop, 23);
            [mut_pop_aa_unique,count_mut_pop_aa,indices_unique] = unique_seqs(mut_pop_aa);
            mut_pop_unique = mut_pop(indices_unique,:);
            mut_pop_bin_unique = mut_pop_bin(indices_unique,:);
            
            %         mut_pop_bin_unique = zeros(size(mut_pop_aa_unique,1),length(H));
            energy_mut_pop_bin_unique = zeros(size(mut_pop_aa_unique,1),1);
            %     energy_mut_pop_bin = [];
            energy_mut_pop_bin = zeros(N,1);
            
            if kk>imm_pres_gen
                bb = b;
            else
                bb = 0;
            end
            
            count_mut_pop = 0;
            count_loop = 1;
            for kkk = 1:size(mut_pop_aa_unique,1)
                %             mut_pop_bin_unique(kkk,:) = convertAAseq2Bin_ahmed(mutant_order,bin_matrix,mut_pop_aa_unique(kkk,:));
                
                %imm_pressure apply or not
                if sum(xor(mut_pop_bin_unique(kkk,imm_pres_site_bin),zeros(1,length(imm_pres_site_bin))),2)==0
                    energy_mut_pop_bin_unique(kkk) = mut_pop_bin_unique(kkk,:)*triu(-H)*mut_pop_bin_unique(kkk,:).' + bb;
                else
                    energy_mut_pop_bin_unique(kkk) = mut_pop_bin_unique(kkk,:)*triu(-H)*mut_pop_bin_unique(kkk,:).';
                    count_mut_pop = count_mut_pop+count_mut_pop_aa(kkk);
                end
                
                %         energy_mut_pop_bin = [energy_mut_pop_bin; repmat(energy_mut_pop_bin_unique(kkk),count_mut_pop_aa(kkk),1)];
                energy_mut_pop_bin(count_loop:count_loop+count_mut_pop_aa(kkk)-1,:) = energy_mut_pop_bin_unique(kkk);
                count_loop = count_loop+count_mut_pop_aa(kkk);
            end
            
            freq_imm_pres_site(kk+1,:) = count_mut_pop/N;
            
            %when freq_imm_pres_site becomes >0 for the first time
            if sum(freq_imm_pres_site(1:kk))==0 && freq_imm_pres_site(kk+1)>0
                escape_mut_appearance = kk+1
            end
            
            Ebar(kk+1) = mean(energy_mut_pop_bin);
            
            %John's metric but does not sum up to 1
            prob_survival_strains = (exp(beta_param*(Ebar(kk+1)-energy_mut_pop_bin_unique)))./...
                (1+exp(beta_param*(Ebar(kk+1)-energy_mut_pop_bin_unique)));
            prob_survival_strains = prob_survival_strains./sum(prob_survival_strains); %for summing to 1
            
            %My modified metric
            %     prob_survival_strains = (exp(beta_param*(Ebar(kk+1)-energy_pop{kk+1})))/sum((exp(beta_param*(Ebar(kk+1)-energy_pop{kk+1}))));
            
            counts_next_gen = mnrnd(N,prob_survival_strains);
            
            %     pop = [];
            %     mut_pop_aa = [];
            %     mut_pop_bin = [];
            pop = char(zeros(N,M));
            mut_pop_aa = char(zeros(N,534));
            mut_pop_bin = char(zeros(N,807));
            count_loop = 1;
            for mm = 1:size(mut_pop_aa_unique,1)
                %         pop = [ pop;  repmat(mut_pop_unique(mm,:),counts_next_gen(mm),1)];
                pop(count_loop:count_loop+counts_next_gen(mm)-1,:) = repmat(mut_pop_unique(mm,:),counts_next_gen(mm),1);
                %         mut_pop_aa = [mut_pop_aa; repmat(mut_pop_aa_unique(mm,:),counts_next_gen(mm),1)];
                mut_pop_aa(count_loop:count_loop+counts_next_gen(mm)-1,:) = repmat(mut_pop_aa_unique(mm,:),counts_next_gen(mm),1);
                %         mut_pop_bin = [mut_pop_bin; repmat(mut_pop_bin_unique(mm,:),counts_next_gen(mm),1)];
                mut_pop_bin(count_loop:count_loop+counts_next_gen(mm)-1,:) = repmat(mut_pop_bin_unique(mm,:),counts_next_gen(mm),1);
                count_loop = count_loop + counts_next_gen(mm);
            end
            
            %         if ismember(kk,(1e3:1e3:1e4))
            %             fprintf('No. of generations = %d\n',kk)
            %         end
            
            if sum(freq_imm_pres_site(kk+1,:) > freq_escape) > 0 %any targeted site > freq_escape
                escape_gen = kk+1
                escape_Ebar = Ebar(kk+1)
                escape_f = freq_imm_pres_site(kk+1,:)
                break;
            end
            
            if kk == Ng %last generation
                escape_gen = kk+1;
                escape_Ebar = Ebar(kk+1);
                escape_f = freq_imm_pres_site(kk+1,:)
            end
            
            %             time_gen(kk) = toc;
            
        end
        
        %     mean(time_gen)
        
        escape_mut_appearance_iter(iter) = escape_mut_appearance
        escape_gen_iter(iter) = escape_gen
        escape_Ebar_iter(iter) = escape_Ebar
        escape_f_iter(iter) = escape_f
    end
    
    escape_gen_iter
    mean(escape_gen_iter)
    escape_f_iter
    
    save(sprintf('results_E1E2/N%d_Ng%d_b%d_beta%.2f_fEscape%.1f_ITER%d_SiteIp%d_Seq%d.mat',...
        N,Ng,b,beta_param,freq_escape,ITER,imm_pres_site,strain_index),...
        'escape_mut_appearance_iter','escape_gen_iter','escape_Ebar_iter','escape_f_iter','indx_seq')
    
end

delete(pobj)