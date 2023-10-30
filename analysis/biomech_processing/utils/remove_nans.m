function [X_out, pid_out, condid_out] = remove_nans(X, pid, condid)
% Helper function to remove rows of nan values in a data matrix and id
% vectors prior to performing PCA analysis.

% Input:
% X = n x m matrix where rows correspond to observations and columns to
% variables. In waveform analysis, this means participant/condition
% indices are by row, and each column is a timepoint in the waveform.
% pid = vector of repeated integers for participant id (1,2,3 etc)
% condid = vector of strings depicting a condition id

% Output:
% X_out, pid_out, condid_out = same data expect any rows with NANs are removed.

    isnanrows = isnan(X);
    pid_out = pid; pid_out(isnanrows(:,1)) = [];
    condid_out = condid; condid_out(isnanrows(:,1))=[];
    X_out = X; X_out(isnanrows(:,1),:) = [];
end