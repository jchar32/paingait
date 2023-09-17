% load the subject-specific trial numbers for each condition collected.

function [subject_info, limb] = get_subject_info(workbookFile)
    subject_info.natural = import_trialnumber(workbookFile, [4,4]);
    subject_info.halfcyc = import_trialnumber(workbookFile, [7,9]);
    subject_info.onecyc = import_trialnumber(workbookFile, [10,12]);
    subject_info.twocyc = import_trialnumber(workbookFile, [13,15]);
    subject_info.threecyc = import_trialnumber(workbookFile, [16,18]);
    subject_info.fourcyc = import_trialnumber(workbookFile, [19,21]);
    subject_info.fivecyc = import_trialnumber(workbookFile, [22,24]);
    subject_info.oneton = import_trialnumber(workbookFile, [25,25]);
    subject_info.threeton = import_trialnumber(workbookFile, [26,26]);
    subject_info.fiveton = import_trialnumber(workbookFile, [27,27]);
    limb = import_trialnumber(workbookFile, [2,2]);
end

function trial_nums = import_trialnumber(workbookFile, dataLines)

    % Set up the Import Options and import the data
    opts = spreadsheetImportOptions("NumVariables", 12);
    
    % Specify sheet and range
    opts.Sheet = "Baseline";
    opts.DataRange = "B" + dataLines(1, 1) + ":M" + dataLines(1, 2);
    
    % Specify column names and types
    opts.VariableNames = ["PAINGAIT01", "PAINGAIT02", "PAINGAIT03", "PAINTGAIT04", "PAINGAIT05", "PAINGAIT06", "PAINGAIT07", "PAINGAIT08", "PAINTGAIT09", "PAINGAIT10", "PAINGAIT11", "PAINGAIT12"];
    opts.VariableTypes = [ "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify variable properties
    % opts = setvaropts(opts, "SubjID", "WhitespaceRule", "preserve");
    % opts = setvaropts(opts, "SubjID", "EmptyFieldRule", "auto");
    
    % Import the data
    trial_nums = readtable(workbookFile, opts, "UseExcel", false);
    
    for idx = 2:size(dataLines, 1)
        opts.DataRange = "A" + dataLines(idx, 1) + ":M" + dataLines(idx, 2);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        trial_nums = [trial_nums; tb]; 
    end

end