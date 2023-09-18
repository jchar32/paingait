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

%% Load participant biomech data

natural

