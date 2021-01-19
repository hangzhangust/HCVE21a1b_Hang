function seq_bin = Binary_Seq(seq,amino_single_combine_array,conserved)
switch nargin
  case 2
    conserved = [];
end
if ~isempty(conserved)
    seq(:,conserved)=[];
end
seq_bin=[];

for i = 1:size(seq,1)
    bin=[];
    for j=1:length(amino_single_combine_array)
        b = char(amino_single_combine_array{j})==seq(i,j);
        if sum(b)==0   
            b(end)=1;
        end
        b = b(2:end);
        b = flip(b);
        bin = [bin b'];
    end
    seq_bin = [seq_bin;bin];
end
end

% E = zeros(size(seq,1),1);
% 
% for i =1:size(seq,1)
% % E(i) = seq_bin(i,:)*J_MPF_BML*seq_bin(i,:)';
% E(i) = seq_bin(i,:)*H';
% end
% % 
% % E = (E-mean(E))/std(E);
