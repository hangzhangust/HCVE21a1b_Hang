
function Gamma = sim_seq(msa)

%Replacing all gaps and ambiguous amino acids in msa by B
msa_int = aa2int(msa);
msa_int(find(msa_int>=21))=21;
msa = int2aa(msa_int);

Code_aminoacid='ACDEFGHIKLMNPQRSTVWYB';

msa_n = zeros(size(msa));
for kk = 1:length(Code_aminoacid)
    msa_n((msa==Code_aminoacid(kk)))=kk;
end

[Nstar,ls]=size(msa_n); 
msa_extended=zeros(Nstar,21*ls);
for kk=1:ls 
    for Num_aminoacid=1:21 
        msa_extended(:,20*(kk-1)+Num_aminoacid)=(msa_n(:,kk)==Num_aminoacid); 
    end
end
Gamma = (msa_extended*msa_extended')/ls; 