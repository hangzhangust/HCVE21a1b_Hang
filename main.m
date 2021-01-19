%%
% *Computational analysis of subtype 1b hepatitis C virus E2 glycoprotein  
% reveals its relatively higher evolutionary flexibility compared to subtype 1a*

% Code for developing the model and re-generating the figures in the paper
%% Setting up paths (of functions and data files required) and necessary parameters

clear;
close all;
clc;

addpath data
addpath functions
addpath 3rd_party_code
addpath figures_scripts

%% Data preparation and model training

% HCV E2 1b sequences were downloaded from GLUE data base
% ((http://hcv-glue.cvr.gla.ac.uk; accessed Oct. 18, 2018)
[header,sequences] = fastaread('AL_1b_E2_aminoAcid_alignment.fasta');

% load patient id and sequence id without patient information
load('NumofPatient_1bE2.mat','No_confidence','NumofPatient')
paient_id  = cell2mat(NumofPatient{1,1});

% find outliers 
outliers = find_outliers(sequences);

% sequence from chimpanzees
chim_id = find(paient_id==0);


% remove outliers, chimpanzee's sequences and sequence without information;
remove_id = unique([chim_id;No_confidence;outliers ]);
sequences(remove_id)=[];
paient_id(remove_id)=[];


% re-weight sequences according to number of sequences from each patient
weight_seq = get_weight_seq(paient_id);

% feed the sequences and weight_seq to the MPF-BML-GUI to train the model
% Code for running MPF-BML is freely available at <https://github.com/ahmedaq/MPF-BML-GUI>. 

%% Fig. 1 Validation of the inferred E2 1b fitness landscape

% (a) plot the energy vs fitness of the inferred model and compare with the
% conservation-only model
% (b) couplings of the model can be used to predict protein contacts
% true contacts were obtained from the crystal structure using Python code
clear;
clc;
run figure1.m

%% Fig. 2 Comparison of the fitness landscapes of E2 1a and 1b

% (a) comparsion of the autocorrelation and neurality associated with each
% landscape. 
% For the calculation of the autocorrelation and neutrality associated with each landscape, 
% check neutrality.m and autocorrelation.m files for details
% (b) peak analysis of each landscape, peaks were generated using peak.m
clear;
clc;
run figure2.m

%% Fig. 3 Validation of E2 1b escape times predicted using the in-host evolutionary model

% (a) comparison of the distribution of escape times associated with the residues 
% at which mutation is known to assist in escape from HmAb
% (b) Comparison of the distribution of escape times associated with the 
% mutations at exposed and buried residues
% escape times were predicting using the evolutionary simulation,
% check WF_E2_supercomp_residue.m for details
clear;
clc;
run figure3.m

%% Fig. 4 Escape time comparison of E2 subtype 1a and 1b

% (a) (c) 
clear;
clc;
run figure4.m
% (d)  
clear;
run figure4_HmAbs.m
% (b) was generated using Pymol, check E21a.pse and E21b.pse for details

%% Supp. Fig. 1 Comparison of the number of parameters required to estimate a fitness landscape for different HCV proteins of subtype 1b
% HCV 1b and 1a E2 seuences (.fasta files, not included in the respository) 
% were downloaded from GLUE data base
% ((http://hcv-glue.cvr.gla.ac.uk; accessed Oct. 20, 2018)
% For the calculation of entropies, check cal_entropy.m
clear;
clc;
run suppfigure1.m

%% Supp. Fig. 2 Statistical validation of the inferred E2 subtype 1b landscape
clear;
clc;
run suppfigure2.m

%% Supp. Fig. 3 Precision of contact predictions vs. the top x pairs obtained using DCA
% the predictions of the DCA were generated using the code dca_Ahmed_W.m, a
% modified version of the original code provided by Morcos et al.
clear;
clc;
run suppfigure3.m

%% Supp. Fig. 4 Robustness tests associated with autocorrelation and neutrality
clear;
clc;
run suppfigure4.m

%% Supp. Fig. 5 & 6 Robustness tests associated with the peaks
% the models were re-inferred following the previous steps using different
% L2 reg. parameters
clear;
clc;
run suppfigure5_and_suppfigure6.m

%% Supp. Fig. 7 Comparison  of  the  residue-wise  entropy  of  E2  subtypes  1a  and  1b
clear;
clc;
run suppfigure7.m

%% Supp. Fig. 8 Comparison  of  predicted  escape  times  of  all  residues  of  E2  1a  and  1b
clear;
clc;
run suppfigure8.m

%% Supp. Fig. 9 Exposed  residues  associated  with  high  escape  time  that  are  commonin  both  subtypes  1a  and  1b
% (a)
clear;
clc;
run suppfigure9.m

% (b) was generated using Pymol, check E21a.pse and E21b.pse for details

%% Supp. Fig. 10 The   B   cell   epitope   coverage   of   E2   for   subtype   1a   and   1b
% epitopes were obtianed from the IEDB database (https://www.iedb.org/; accessed Jan. 7, 2021)
clear;
clc;
run suppfigure10.m

% (b) was generated using Pymol, check E21a.pse and E21b.pse for details

%% Supp. Fig. 11 Designing  a  classifier  to  determine  optimal  escape  time  for  E2  1bbased on the available knowledge of experimentally or clinically identified escapemutations.
clear;
clc;
run suppfigure11.m

% (b) was generated using Pymol, check E21a.pse and E21b.pse for details