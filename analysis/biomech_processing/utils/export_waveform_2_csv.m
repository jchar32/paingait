% export waveform data to csv

load(fullfile("../data", "waveform_data.mat"), "wfrm")

outcomes = ["hip","knee","ankle", "grf"];
comparisons = ["cyca","cycb","cycc","cycavton","cyc","ton"];
allcondnames = fieldnames(wfrm.hip.l.angle);

for o = 1:length(outcomes)
    vars = ["angle","moment","force"];
    for v = 1:2
        if strcmp(outcomes{o}, "grf")
            vname = vars(3);
        else
            vname = vars(v);
        end

        matrix_tbl = table();
        matrix_mean_tbl = table();
        
        for c = 1:length(allcondnames)
            data = wfrm.(outcomes{o}).l.(vname).(allcondnames{c});
            datamean(:,1) = nanmean(data(:,:,1),1);
            datamean(:,2) = nanmean(data(:,:,2),1);
            datamean(:,3) = nanmean(data(:,:,3),1);

            data = [reshape(data(:,:,1)', 12*101,1), ...
                reshape(data(:,:,2)', 12*101,1), ...
                reshape(data(:,:,3)', 12*101,1)];
            
            matrix = [];
            matrix_mean = [];
            for a = 1:3 % axis
                matrix = [matrix;[[data(:,a), repmat([0:1:100]',12,1)], ...
                    string(repelem([1:1:12]',101,1)), ...
                    string(repmat([allcondnames{c}],12*101,1)), ...
                    string(repmat("l", 12*101,1)), ...
                    string(repmat(num2str(a),12*101,1))]];
                matrix_mean = [matrix_mean;[[datamean(:,a), [0:1:100]'], ...
                    string(repelem([1]',101,1)), ...
                    string(repmat([allcondnames{c}],101,1)), ...
                    string(repmat("l", 101,1)), ...
                    string(repmat(num2str(a),101,1))]];
            end

            if c == 1
                matrix_tbl = array2table(matrix, variablenames=["val", "perc","pid","cond","limb", "axis"]);
                matrix_tbl = convertvars(matrix_tbl,"val", "double");
                matrix_tbl = convertvars(matrix_tbl,["pid","perc","cond","limb" "axis"], "categorical");

                matrix_mean_tbl = array2table(matrix_mean, variablenames=["val", "perc","pid","cond","limb", "axis"]);
                matrix_mean_tbl = convertvars(matrix_mean_tbl,"val", "double");
                matrix_mean_tbl = convertvars(matrix_mean_tbl,["pid","perc","cond","limb" "axis"], "categorical");
            else
                matrix_tbl_temp = array2table(matrix, variablenames=["val", "perc","pid","cond","limb", "axis"]);
                matrix_tbl_temp = convertvars(matrix_tbl_temp,"val", "double");
                matrix_tbl_temp = convertvars(matrix_tbl_temp,["pid","perc","cond","limb" "axis"], "categorical");
                matrix_tbl = [matrix_tbl; ...
                    matrix_tbl_temp];                

                matrix_mean_tbl_temp = array2table(matrix_mean, variablenames=["val", "perc","pid","cond","limb", "axis"]);
                matrix_mean_tbl_temp = convertvars(matrix_mean_tbl_temp,"val", "double");
                matrix_mean_tbl_temp = convertvars(matrix_mean_tbl_temp,["pid","perc","cond","limb" "axis"], "categorical");
                matrix_mean_tbl = [matrix_mean_tbl; ...
                    matrix_mean_tbl_temp];

            end
        end
        
        writetable(matrix_tbl, fullfile("../data/waveforms",[(outcomes{o}) + "_" + vname + ".csv"] ));
        writetable(matrix_mean_tbl, fullfile("../data/waveforms",["mean" + (outcomes{o}) + "_" + vname + ".csv"] ))

    end
end


