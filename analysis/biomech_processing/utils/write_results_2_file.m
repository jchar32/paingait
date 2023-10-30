function write_results_2_file(S, filename)
% write a results table to a csv file.

% Input:
% S = structure of PCA analysis results (arising from pca_svd function).
% filename = the file name to be written to.

% Output:
% None

    varnames = [repelem("PC",S.n_pcs2keep) + string(1:1:S.n_pcs2keep)];
    results_table = array2table(S.testablepcscores, VariableNames=varnames);
    results_table = addvars(results_table, S.pid, S.condid, NewVariableNames=["pid","condid"]);  
    writetable(results_table, fullfile("../data/pca/", filename));
end