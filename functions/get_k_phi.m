function k = get_k_phi(Profile,phi_opt)
 k = zeros(1,size(Profile,2));
for i =1:size(Profile,2)
   fre = sort(Profile(:,i),'descend');
   fre = fre(fre>0);
   entropy = -sum(fre.*log(fre));
   if length(fre)>1
       entropy_ratio = zeros(length(fre)-1,1);
       for ki=1:length(fre)-1
        fbar = sum(fre(ki+1:end)); % frequency of all amino acids
        entropy_ki = -sum(fre(1:ki).*log(fre(1:ki))) - fbar*log(fbar);
        entropy_ratio(ki) = entropy_ki/entropy;
       end
       ki_phi = find(entropy_ratio>=phi_opt,1);
       if length(ki_phi>0)
       k(i) = ki_phi;
       else
           k(i) = length(fre)-1;
       end
   end
   
end