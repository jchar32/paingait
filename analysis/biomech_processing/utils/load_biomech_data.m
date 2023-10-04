% loads in a participants dataset for all conditions in study
function [data_out] = load_biomech_data(p, data_root, subject_info)

% convert participant number to 0X format
partic = num2str(p);
if p<10; partic = "0"+num2str(p); end

% load eahc .mat file and place contents into a struct with a field name based on condition
% --cyc conditions had 3 trials which get labled a,b,c
condition_names = fieldnames(subject_info);
for t = 1:length(condition_names) % 23 total trials to compile
   
    % load .mat file for each condition
    if size(subject_info.(condition_names{t}),1) > 1
        
        increment = ["a","b","c"]; % suffix to add to the 1-3 trials in cyc condition
        for subt = 1:size(subject_info.(condition_names{t}),1)
            
            t_num = subject_info.(condition_names{t}).("PAINGAIT" + partic)(subt);
            % check if this condition wasnt completed
            if isnan(t_num)
                data_out.(condition_name_updated) = [];
                continue;
            end
            
            % Cortex 8 adds a "0" ahead of single digit trial numbers, this accomodates it. Cortex 8 was used from P10 onward I believe.
            if exist(fullfile(data_root, "PAINGAIT" + partic, "Trial" + num2str(t_num) + ".mat"),"file") == 0
                path2file = fullfile(data_root, "PAINGAIT" + partic, "Trial0" + num2str(t_num) + ".mat");
            else
                path2file = fullfile(data_root, "PAINGAIT" + partic, "Trial" + num2str(t_num) + ".mat");
            end
            
            condition_name_updated = condition_names{t} + increment(subt);
            
            [data_out.(condition_name_updated)] = load(path2file);
        end
    else
        
        % check if this condition wasnt completed
        if any(isnan(subject_info.(condition_names{t}).("PAINGAIT" + partic)(1)))
            data_out.(condition_names{t}) = [];
            continue;
        end

        path2file = fullfile(data_root, "PAINGAIT" + partic, "Trial" + num2str(subject_info.(condition_names{t}).("PAINGAIT" + partic)(1)) + ".mat");
        
        [data_out.(condition_names{t})] = load(path2file);
    end
     
end

end
