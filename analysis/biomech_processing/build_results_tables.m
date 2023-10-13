function [outtbl] = build_results_tables(discrete_data, p)
outtbl=struct();
heads = fieldnames(discrete_data);
subcats = fieldnames(discrete_data.(heads{1}));

for c = 1:size(fieldnames(discrete_data),1)

    for subcat = 1:size(subcats,1)
        if subcat ~= 1
            % Right
            outcomes_time = discrete_data.(heads{c}).(subcats{subcat}).r.means.time;
            outcomes_time.limb="r";
            outcomes_time.cond = string(heads{c});
            outcomes_time.pid = p;
            outtbl.(subcats{subcat}).time_outcomes_r(c,:) = struct2table(outcomes_time);

            outcomes = discrete_data.(heads{c}).(subcats{subcat}).r.means;
            outcomes_notime = rmfield(outcomes,'time');
            outcomes_notime.limb="r";
            outcomes_notime.cond = string(heads{c});
            outcomes_notime.pid = p;
            outtbl.(subcats{subcat}).biomech_outcomes_r(c,:) = struct2table(outcomes_notime);

            % Left
            outcomes_time = discrete_data.(heads{c}).(subcats{subcat}).l.means.time;
            outcomes_time.limb="l";
            outcomes_time.cond = string(heads{c});
            outcomes_time.pid = p;
            outtbl.(subcats{subcat}).time_outcomes_l(c,:) = struct2table(outcomes_time);

            outcomes = discrete_data.(heads{c}).(subcats{subcat}).l.means;
            outcomes_notime = rmfield(outcomes,'time');
            outcomes_notime.limb="l";
            outcomes_notime.cond = string(heads{c});
            outcomes_notime.pid = p;
            outtbl.(subcats{subcat}).biomech_outcomes_l(c,:) = struct2table(outcomes_notime);

        else
            % Right
            outcomes = discrete_data.(heads{c}).(subcats{subcat}).r.means;
            outcomes.limb="r";
            outcomes.cond = string(heads{c});
            outcomes.pid = p;
            outtbl.(subcats{subcat}).biomech_outcomes_r(c,:) = struct2table(outcomes);

            % Left
            outcomes = discrete_data.(heads{c}).(subcats{subcat}).l.means;
            outcomes.limb="l";
            outcomes.cond = string(heads{c});
            outcomes.pid=p;
            outtbl.(subcats{subcat}).biomech_outcomes_l(c,:) = struct2table(outcomes);
        end
    end
end

end