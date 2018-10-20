function [Z_mat,par_child_mat ] = Z_mat_sim(S,C)


Z_mat = zeros(S,C);

S_perm = randperm(S);

p = (1/C)*ones(1,C);

r = mnrnd(S,p);

indic = length(find(r == 0)); %% i.e non has 0 mutation in it

if indic ~= 0 
    ind_zer = find(r == 0);
    
    r(ind_zer) = 1; % assign 1mut each to the ones without mut
    
    [val,ind_max] = max(r);
    
    r(ind_max) = val - length(ind_zer);  % subtract the added mut from the max  
end


    rr = cumsum(r);
    rrr = [0 rr];
    

    for c = 1:C

        beg_c = rrr(c) + 1;
        end_c = rrr(c+1);

        for h = beg_c:end_c
            
            Z_mat(S_perm(h),c) = 1;
        end   
    end


%Z_mat

C_perm = randperm(C);

for c = 1:C
    
    
    if c == 1  % first clone 
       c_per = C_perm(1); 
    else
        
        new_child = C_perm(c);
        %%% possible parents
        possible_parents = C_perm(1:c-1);
        
        ind_chosen_par = randperm(length(possible_parents),1);
        
        chosen_par = possible_parents(ind_chosen_par);
        
        %%%% all the muts of the parent goes to the child
        
        chosen_par_col = Z_mat(:,chosen_par); 
        new_clone_col = Z_mat(:,new_child);
        
        ind_par_mut = find(chosen_par_col == 1);
        
        new_clone_col(ind_par_mut) = 1;
        
        Z_mat(:,new_child) = new_clone_col;
        
        par_child_mat(c-1,:) = [chosen_par new_child];
        
    end
   
end


end

