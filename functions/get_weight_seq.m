function weight_seq = get_weight_seq(patient)
C = tabulate(patient);
weight_seq = zeros(length(patient),1);

for i =1:length(weight_seq)
    for j =1:size(C,1)
        if(patient(i)==C(j,1))
            weight_seq(i) = 1/C(j,2);
        end
    end
end
end