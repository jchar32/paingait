


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