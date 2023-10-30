function [pid, condid] = get_data_ids(X, allconds, comps)
% based on the comparisons of interest, get the conditions needed for that
% comparisons.
% A Helper function for compiling data for PCA
    if strcmp(comps,"cycavton")
        conds = allconds(contains(allconds,["onecyca","threecyca","fivecyca","ton"]));
    else
        conds = allconds(contains(allconds, comps));
    end

    condid = repelem(conds,size(X,1)/length(conds));
    pid = repmat([1:1:12]',[size(X,1)/12,1]);
end