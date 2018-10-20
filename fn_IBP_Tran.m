function [Bin_Mat_Loni] = fn_IBP_Tran(Bin_Mat_Lana,alpha,time_step)


if time_step == 1
    
    K_i = poissrnd(alpha); 
    
    Bin_Mat_Loni = ones(1,K_i);
    
    
    
    
else
   
    non_zero_col = zeros(1,size(Bin_Mat_Lana,2));
    
    for j = 1:size(Bin_Mat_Lana,2)
        col_j = Bin_Mat_Lana(:,j);
        
        if sum(col_j) > 0
           non_zero_col(j) = 1; 
        end
    end
    
    
    K_plus = sum(non_zero_col);
    
    new_row_part_1 = zeros(1,K_plus);
    
    for k = 1:K_plus
        kkk = Bin_Mat_Lana(:,k);
        m_ik = sum(kkk);
        
        p_success = m_ik/time_step;
        
        new_row_part_1(k) = binornd(1,p_success); 
        
    end
    
    
    K_i_new = poissrnd(alpha/time_step);
    
            if K_i_new == 0
                
                Bin_Mat_Loni = [Bin_Mat_Lana;new_row_part_1];

            else
                new_row_part_2 = ones(1,K_i_new);
                new_row = [new_row_part_1 new_row_part_2];


                padding_zeros = zeros(size(Bin_Mat_Lana,1),K_i_new);

                Bin_Mat_Loni = [Bin_Mat_Lana  padding_zeros; new_row];

            end
    
end

end

