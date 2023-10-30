function [out] = pca_svd(X)
% calculate principal components, scores, and variance explained using the
% singular value decomposition approach.

% Input:    
% X: nxm numeric matrix where n=observations/participants/participant-conditions & m=variables/timepoints/samplepoints
% Output:
% out: struct containing the PCs, PCscores, eigenvalues, variance explained
%       and the number of pcs to keep based on a 90% threshold as used in common
%       biomechancis literature (e.g., Hatfield et al 2015 DOI:
%       10.1002/acr.22564)

    Xmean_cent = X - mean(X).*ones(size(X,1),1);
    [U,S,Vt] = svd(Xmean_cent, 'econ');
    pcs = Vt;
    eigvals = (diag(S).^2) ./ (size(Xmean_cent,1)-1);
    pcscores = Xmean_cent*pcs;
    varexp = eigvals / sum(eigvals);
    n_pcs_2keep = find(cumsum(varexp) > 0.9, 1,'first');
    
    testablepcscores = pcscores(:,1:n_pcs_2keep);
    out.pcs = pcs; out.pcscores=pcscores; out.eigvals=eigvals; out.varexp = varexp;
    out.n_pcs2keep=n_pcs_2keep; out.testablepcscores = testablepcscores;
end

