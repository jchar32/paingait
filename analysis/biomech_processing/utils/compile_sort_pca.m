function [pc_percentiles_vec] = compile_sort_pca(res, X)
    pc_percentiles_vec = [];
    varnames = [repelem("PC",res.n_pcs2keep) + string(1:1:res.n_pcs2keep)];
    pc = array2table(res.testablepcscores, VariableNames=varnames);
    pc = addvars(pc, res.pid, res.condid, NewVariableNames=["pid","condid"]);
    
    pc1_sorted = sortrows(pc,"PC1","ascend");
    [pc_percentiles_vec] = index_5th_95th_percentile(pc1_sorted, pc, res, pc_percentiles_vec, X);
%     
%     pc_percentiles = [round(size(pc,1)*0.05), round(size(pc,1)*0.95)];
%     pc_percentiles_idx = pc1_sorted(pc_percentiles,size(pc,2)-1:size(pc,2));
%     
%     rowids1 = [find(res.pid == pc_percentiles_idx.pid(1) & strcmp(res.condid, pc_percentiles_idx.condid(1))), ...
%             find(res.pid == pc_percentiles_idx.pid(2) & strcmp(res.condid, pc_percentiles_idx.condid(2)))];
%     pc_percentiles_vec(:,:) = [pc_percentiles_vec;X(rowids1,:)];
    
    if length(varnames) > 1
        pc2_sorted = sortrows(pc,"PC2","ascend");
        [pc_percentiles_vec] = index_5th_95th_percentile(pc2_sorted, pc, res, pc_percentiles_vec, X);
%         
%         pc_percentiles_idx = pc1_sorted(pc_percentiles,size(pc,2)-1:size(pc,2));
%         rowids2 = [find(res.pid == pc_percentiles_idx.pid(1) & strcmp(res.condid, pc_percentiles_idx.condid(1))), ...
%                 find(res.pid == pc_percentiles_idx.pid(2) & strcmp(res.condid, pc_percentiles_idx.condid(2)))];
%         pc_percentiles_vec(:,:) = [pc_percentiles_vec;X(rowids2,:)];
    end


end

function [pc_percentiles_vec] = index_5th_95th_percentile(pcsorted, pc, res, pc_percentiles_vec, X)
    pc_percentiles = [round(size(pc,1)*0.05), round(size(pc,1)*0.95)];
    pc_percentiles_idx = pcsorted(pc_percentiles, size(pc,2)-1:size(pc,2));
    
    rowids = [find(res.pid == pc_percentiles_idx.pid(1) & strcmp(res.condid, pc_percentiles_idx.condid(1))), ...
            find(res.pid == pc_percentiles_idx.pid(2) & strcmp(res.condid, pc_percentiles_idx.condid(2)))];
    pc_percentiles_vec = [pc_percentiles_vec; X(rowids,:)];
end