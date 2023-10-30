function [pca_data] = gather_wfrm_comparisons(wfrm, outcomes, comparisons)
% Compile waveform data into matrices of various comparisons for PCA analysis

% Input: 
% wfrm: struct containing the row-wise waveform data in an heierarch of
%   outcome (e.g., hip) --
%       limb (e.g., r or l) --
%         variable (e.g., angle) --
%           condition (e.g., natural, onecyca) --
%                 data = [participant x samplepoint x axis] matrix (e.g.,
%                 12 x 100 x 3)
% outcomes = list of strings corresponding to the outcomes to be compiled
% comparisons = list of strings corresponding to the condition comparisons
% to be compiled
% Output:
%     pca_data: structure with similar hierarchy to wfrm but conditions
%     have been gathered and vertically stacked. The resulting matrix holds
%     (row wise) each participant/condtion combination for a given set of
%     comparisons.
    for o = 1:length(outcomes)
    
        for c = 1:length(comparisons)
            comps = comparisons{c};
            if strcmp(outcomes{o}, "grf")
                [id] = get_condition_ids(comps, wfrm.(outcomes{o}).l.force);
                temp = struct2cell(wfrm.(outcomes{o}).l.force);
                temp=temp(id);
                pca_data.(outcomes{o}).l.(comps).force = vertcat(temp{:,:,:});
            
                temp = struct2cell(wfrm.(outcomes{o}).r.force);
                temp=temp(id);
                pca_data.(outcomes{o}).r.(comps).force = vertcat(temp{:,:,:});
            else
                [id] = get_condition_ids(comps, wfrm.(outcomes{o}).l.angle);
                temp = struct2cell(wfrm.(outcomes{o}).l.angle);
                temp=temp(id);
                pca_data.(outcomes{o}).l.(comps).angle = vertcat(temp{:,:,:});
            
                temp = struct2cell(wfrm.(outcomes{o}).r.angle);
                temp=temp(id);
                pca_data.(outcomes{o}).r.(comps).angle = vertcat(temp{:,:,:});
        
                temp = struct2cell(wfrm.(outcomes{o}).l.moment);
                temp=temp(id);
                pca_data.(outcomes{o}).l.(comps).moment = vertcat(temp{:,:,:});
            
                temp = struct2cell(wfrm.(outcomes{o}).r.moment);
                temp=temp(id);
                pca_data.(outcomes{o}).r.(comps).moment = vertcat(temp{:,:,:});
            end
            
        end
    
    end

end

function [ids, condids] = get_condition_ids(comps, condstruct)

    if strcmp(comps, "cyca")
        ids = contains(fieldnames(condstruct), ["natural", "cyca"]);
    elseif strcmp(comps, "cycb")
        ids = contains(fieldnames(condstruct), ["natural","cycb"]);
    elseif strcmp(comps, "cycc")
        ids = contains(fieldnames(condstruct), ["natural","cycc"]);
    elseif strcmp(comps, "cycavton")
        ids = contains(fieldnames(condstruct), ["natural","onecyca","threecyca","fivecyca","ton"]);
    elseif strcmp(comps, "cyc")
        ids = contains(fieldnames(condstruct), ["natural","cyc"]);
    elseif strcmp(comps, "ton")
        ids = contains(fieldnames(condstruct), ["natural","ton"]);
    end
end
