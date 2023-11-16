function [pc_percentiles_vec] = compile_sort_pca(res, X)
    pc_percentiles_vec = [];
    varnames = [repelem("PC",res.n_pcs2keep) + string(1:1:res.n_pcs2keep)];
    pc = array2table(res.testablepcscores, VariableNames=varnames);
    pc = addvars(pc, res.pid, res.condid, NewVariableNames=["pid","condid"]);
    
    pc1_sorted = sortrows(pc,"PC1","ascend");
    [pc_percentiles_vec] = index_5th_95th_percentile(pc1_sorted, pc, res, pc_percentiles_vec, X);
    
    if length(varnames) > 1
        pc2_sorted = sortrows(pc,"PC2","ascend");
        [pc_percentiles_vec] = index_5th_95th_percentile(pc2_sorted, pc, res, pc_percentiles_vec, X);

    end
end

function [pc_percentiles_vec] = index_5th_95th_percentile(pcsorted, pc, res, pc_percentiles_vec, X)
    pc_percentiles = [round(size(pc,1)*0.05), round(size(pc,1)*0.95)];
    pc_percentiles_idx = pcsorted(pc_percentiles, size(pc,2)-1:size(pc,2));
    
    lower_5th = 1:round(size(pc,1)*0.05);
    upper_5th = round(size(pc,1)*0.95:size(pc,1));
    pc_perc_idx_lower = pcsorted(lower_5th, size(pc,2)-1:size(pc,2));
    pc_perc_idx_upper = pcsorted(upper_5th, size(pc,2)-1:size(pc,2));
    
    for i = 1:size(pc_perc_idx_lower,1)
        row_id_lower(i) = find(res.pid == pc_perc_idx_lower.pid(i) & strcmp(res.condid, pc_perc_idx_lower.condid(i)));
    end
    for i = 1:size(pc_perc_idx_lower,1)
        row_id_upper(i) = find(res.pid == pc_perc_idx_upper.pid(i) & strcmp(res.condid, pc_perc_idx_upper.condid(i)));
    end
    
    mean_wfrm_lower = nanmean(X(row_id_lower,:));
    mean_wfrm_upper = nanmean(X(row_id_upper,:));

    pc_percentiles_vec = [pc_percentiles_vec; mean_wfrm_lower; mean_wfrm_upper];
   
end