% Processing for PAINGAIT study

% 1. ingest trial numbers for participants that correspond to each walking trial and condition
% 2. identify gait cycles to use for analysis
% 3. calculate gait outcomes
% 4. compile waveforms and time normalize


%% Setup


data_direc_root = "U:\Projects\Experimental pain gait assessment\";

subjec_info = get_subject_info(fullfile(data_direc_root,"Processed Data","Experimental pain gait assessment subject trial info.xlsx"));



