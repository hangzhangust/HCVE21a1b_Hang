%% Peak Analysis
load('Only_confidence_best.mat')
Peaks = [];
tic;
num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
for i =1:size(msa_aa,1)
    
    x = msa_bin(i,:);

    x_E=x*J_MPF_BML*x';

    while 1
        lowest_E=x_E;
        for j =1:length(num_mutants_combine_array)
            Start_site = num_mutants_combine_array_acc_all(j)+1;
            End_site   = num_mutants_combine_array_acc_all(j+1);
            one = find(x(Start_site:End_site)>0);
            if(one)
                zero_E = x_E-2*x*J_MPF_BML(:,Start_site+one-1)+J_MPF_BML(Start_site+one-1,Start_site+one-1);

            else
                zero_E = x_E;
            end
            low_E=zero_E;
            pos_start=Start_site;
            pos_end=End_site;
            pos_change=-1;
            for k = Start_site:End_site
                tmp_x=x;
                tmp_x(Start_site:End_site)=0;
                tmp_x(k)=1;
                E = zero_E+2*tmp_x*J_MPF_BML(:,k)-J_MPF_BML(k,k);
                if E<low_E
                    low_E = E;

                    pos_change=k;
                end
            end

            if low_E<lowest_E
                lowest_E = low_E;
                pos_start_all=pos_start;
                pos_end_all=pos_end;
                pos_change_all=pos_change;
            end
        end

        if abs(x_E-lowest_E)>1e-7
            x_E=lowest_E;
            if pos_change_all>0
               x(pos_start_all:pos_end_all)=0;
               x(pos_change_all)=1;
               
            else
               x(pos_start_all:pos_end_all)=0;
            end
        else
            Peaks=[Peaks;x];
            break;
        end
        
    end
    

end
