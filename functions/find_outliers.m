function outliers = find_outliers(sequences)
%%

for kk = 1:length(sequences)
    msa(kk,:) = sequences{kk};
end
S = sim_seq(msa);

[eigvec,~] = eig_sort(S);

%% Outliers

o1 = find(isoutlier(eigvec(:,1)));
o2 = find(isoutlier(eigvec(:,2)));
% o3 = find(isoutlier(eigvec(:,3)));
outliers = unique([o1;o2]);





