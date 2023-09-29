% Processing for PAINGAIT study

% 1. ingest trial numbers for participants that correspond to each walking trial and condition
% 2. identify gait cycles to use for analysis
% 3. calculate gait outcomes
% 4. compile waveforms and time normalize


%% Setup
addpath("./utils")

data_direc_root = "U:\Projects\Experimental pain gait assessment\";
data_root = data_direc_root + "\Processed Data\Motion Analysis Data\";
[subject_info, limb] = get_subject_info(fullfile(data_direc_root,"Processed Data","Experimental pain gait assessment subject trial info.xlsx"));

sample_rate.mocap= 100;
sample_rate.analog= 2000;

%% Load participant biomech data and process

for p = 1:size(subject_info.natural,2)
    all_data = load_biomech_data(p, data_root,subject_info);
    [gait_events] = process_gait_events(all_data);
    [discrete_data] = calculate_discrete_outcomes(all_data, gait_events, sample_rate);
    [waveforms] = compile_waveforms(all_data, gait_events, sample_rate);
    visualize_gait_waveforms(waveforms)
end

