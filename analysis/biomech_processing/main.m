% Processing for PAINGAIT study

% 1. ingest trial numbers for participants that correspond to each walking trial and condition
% 2. identify gait cycles to use for analysis
% 3. calculate gait outcomes
% 4. compile waveforms and time normalize
% 5. gather data into results table for analysis

%% Setup
addpath("./utils")

data_direc_root = "U:\Projects\Experimental pain gait assessment\";
data_root = data_direc_root + "\Processed Data\Motion Analysis Data\";
[subject_info, limb] = get_subject_info(fullfile(data_direc_root,"Processed Data","Experimental pain gait assessment subject trial info.xlsx"));

sample_rate.mocap= 100;
sample_rate.analog= 2000;

%% Load participant biomech data and process

for p = 12:size(subject_info.natural,2)
    all_data = load_biomech_data(p, data_root,subject_info);
    [gait_events] = process_gait_events(all_data);
    [discrete_data] = calculate_discrete_outcomes(all_data, gait_events, sample_rate);
    [waveforms] = compile_waveforms(all_data, gait_events, sample_rate);
    visualize_gait_waveforms(waveforms)
    [results_table] = build_results_tables(discrete_data, p);
    
    uiwait(msgbox("Look through participant waveforms, note any anomalies. Click OK when done to save data")); close all;
    save(fullfile("../data", ["P" + num2str(p) + "_data.mat"]), "gait_events","discrete_data","waveforms", "results_table")
    disp("Done: " + num2str(p))
end

%% Compile participant data into summary tables for analysis

temporal = table();
knee = table();
hip = table();
ankle = table();
grf = table();
thigh = table();
shank = table();
foot = table();

for p=1:size(subject_info.natural,2)
    load(fullfile("../data", ["P" + num2str(p) + "_data.mat"]), "waveforms", "results_table");
    temporal = [temporal; results_table.temporal.biomech_outcomes_l;results_table.temporal.biomech_outcomes_r];
    knee = [knee; results_table.knee.biomech_outcomes_l;results_table.knee.biomech_outcomes_r];
    hip = [hip; results_table.hip.biomech_outcomes_l;results_table.hip.biomech_outcomes_r];
    ankle = [ankle; results_table.ankle.biomech_outcomes_l;results_table.ankle.biomech_outcomes_r];
    grf = [grf; results_table.grf.biomech_outcomes_l;results_table.grf.biomech_outcomes_r];
    thigh = [thigh; results_table.thigh.biomech_outcomes_l;results_table.thigh.biomech_outcomes_r];
    shank = [shank; results_table.shank.biomech_outcomes_l;results_table.shank.biomech_outcomes_r];
    foot = [foot; results_table.foot.biomech_outcomes_l;results_table.foot.biomech_outcomes_r];
end

writetable(temporal,fullfile("../data","temporal.csv"))
writetable(knee,fullfile("../data","knee.csv"))
writetable(hip,fullfile("../data","hip.csv"))
writetable(ankle,fullfile("../data","ankle.csv"))
writetable(grf,fullfile("../data","grf.csv"))
writetable(thigh,fullfile("../data","thigh.csv"))
writetable(shank,fullfile("../data","shank.csv"))
writetable(foot,fullfile("../data","foot.csv"))

%% Compile waveforms into matrices for analysis

for p=1:size(subject_info.natural,2)
    load(fullfile("../data", ["P" + num2str(p) + "_data.mat"]), "waveforms");
    condition_names = fieldnames(waveforms);
    for c = 1:size(condition_names)
        try waveforms.(condition_names{c});
        catch
            continue; % skip condition as it was not collected
        end
        outcome_names = fieldnames(waveforms.(condition_names{c}));
        for o = 1:size(outcome_names,1)
            limbs =["l","r"];
            for limb = 1:2
                varnames = fieldnames(waveforms.(condition_names{c}).(outcome_names{o}).(limbs{limb}));
                for v = 1:size(varnames)
                    if strcmp(varnames{v}, "force") | strcmp(varnames{v}, "angle")
                        gait_comp = "cycle_nd_mean";
                    else
                        gait_comp = "stance_nd_mean";
                    end

                    if p ==1 % initialize array with nans so they don't default to zeros
                        wfrm.(outcome_names{o}).(limbs{limb}).(varnames{v}).(condition_names{c})= nan(12,100,3);
                    end

                    wfrm.(outcome_names{o}).(limbs{limb}).(varnames{v}).(condition_names{c})(p,1:100,1:3) = ...
                        waveforms.(condition_names{c}).(outcome_names{o}).(limbs{limb}).(varnames{v}).(gait_comp);
                end
            end
        end
    end
    clear waveforms
end

save(fullfile("../data", "waveform_data.mat"), "wfrm")

%% Compile data and PCA

% compile various combinations of conditions for PCA comparisons
% 1. all
% 2. start
% 3. mid
% 4. end
% 5. cyc start 1,2,3 vs tonic
% 6. all cyc

% Outcomes
% hip x,y,z angle, moment
% knee x,y,z angle, moment
% ankle x,y,z angle moment
% grf x,y,z force

load(fullfile("../data", "waveform_data.mat"), "wfrm")

outcomes = ["hip","knee","ankle", "grf"];
comparisons = ["cyca","cycb","cycc","cycavton","cyc","ton"];
allcondnames = fieldnames(wfrm.hip.l.angle);

% for each comparisons of interest, gather together participant/condition
% data into single matricies for analysis
[pca_data] = gather_wfrm_comparisons(wfrm, outcomes, comparisons);

% Really ugly nested loops to step down into hierearchy, get specific data,
% run pca, store results. There is probably a more elegant way but this is
% fine...
for o = 1:length(outcomes)
    for c = 1:length(comparisons)
        for lb = 1:2
            if lb == 1; limb = "r"; else; limb = "l"; end
            
            for axis = 1:3
                if strcmp(outcomes{o},"grf")
                    % get data and make id columns for analysis
                    X=pca_data.(outcomes{o}).(limb).(comparisons{c}).force(:,:,axis);
                    [pid, condid] = get_data_ids(X, allcondnames, comparisons{c});
                    [X, pid, condid] = remove_nans(X, pid, condid);
                    
                    % perform pca on waveforms
                    [S] = pca_svd(X);
                    S.pid=pid; S.condid=condid; S.X = X;
                    pca_results.(outcomes{o}).(limb).(comparisons{c}).force{axis} = S;

                    % generate file name
                    filename = strjoin([outcomes{o}, limb, comparisons{c}, "force", num2str(axis) + ".csv"],"_");
                    write_results_2_file(S, filename)
                else
                    % ANGLE
                    % get data and make id columns for analysis
                    X=pca_data.(outcomes{o}).(limb).(comparisons{c}).angle(:,:,axis);
                    [pid, condid] = get_data_ids(X, allcondnames, comparisons{c});
                    [X, pid, condid] = remove_nans(X, pid, condid);
                    
                    % perform pca on waveforms
                    [S] = pca_svd(X);
                    S.pid=pid; S.condid=condid; S.X = X;
                    pca_results.(outcomes{o}).(limb).(comparisons{c}).angle{axis} = S;

                    % generate file name
                    filename = strjoin([outcomes{o}, limb, comparisons{c}, "angle", num2str(axis) + ".csv"],"_");
                    write_results_2_file(S, filename)
                    
                    % MOMENT
                    % get data and make id columns for analysis
                    X=pca_data.(outcomes{o}).(limb).(comparisons{c}).moment(:,:,axis);
                    [pid, condid] = get_data_ids(X, allcondnames, comparisons{c});
                    [X, pid, condid] = remove_nans(X, pid, condid);
                    
                    % perform pca on waveforms
                    [S] = pca_svd(X);
                    S.pid=pid; S.condid=condid; S.X = X;
                    pca_results.(outcomes{o}).(limb).(comparisons{c}).moment{axis} = S;

                    % generate file name
                    filename = strjoin([outcomes{o}, limb, comparisons{c}, "moment", num2str(axis) + ".csv"],"_");
                    write_results_2_file(S, filename)
                    
                end % If statement
                
            end % axis

        end % limb
    end % comparison
end % outcome

save(fullfile("../data", "pca_analysis.mat"), "pca_data","pca_results","comparisons","outcomes")

%% Creating figures for pca analysis

load(fullfile("../data", "pca_analysis.mat"), "pca_data","pca_results","comparisons","outcomes")

% plot pca results
varname = ["angle","moment","force"];
for o = 1:length(outcomes)
    
    for c = 1:length(comparisons)
        
        for v = 1:3
            if strcmp(outcomes{o}, "grf")
                vname = varname(3);
            else
                if v==3; continue;end
                vname = varname(v);
            end

            t=tiledlayout(2,3);
            t.TileSpacing="compact";
            t.Padding="compact";
            title(t,strjoin([outcomes{o}, comparisons{c}, vname]," "));
            for axis = 1:3
                % extract data
                res=pca_results.(outcomes{o}).l.(comparisons{c}).(vname){axis};
                X = res.X;
                Xmean = nanmean(X);
                
                % calc means for each condition in dataset
                comps = unique(res.condid,"stable");
                for i = 1:length(comps)
                    condmeans(:,i) = nanmean(X(strcmp(res.condid, comps{i}),:),1)';
                end
                condmeanstbl = array2table(condmeans,VariableNames=comps);
                [pc_percentiles_vec] = compile_sort_pca(res, X);

                % plot the compiled pca results waveforms
                nexttile(axis); hold on;
                pca_plot(X, Xmean, pc_percentiles_vec, axis)
                
                nexttile(axis+3); hold on;
                mean_wfrm_plot(Xmean, condmeanstbl, axis);
                clear condmeans condmeanstbl pc_percentiles_vec
            end
            
        filename = strjoin([outcomes{o}, "l", comparisons{c}, vname + ".fig"],"_");
        savefig(fullfile("../data/pca/figs", filename ));
        close all
        end
        
    end
    
end












