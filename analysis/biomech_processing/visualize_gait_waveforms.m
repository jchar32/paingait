function visualize_gait_waveforms(data)
% plot several key waveforms to quality checking purposes.

conditions = fieldnames(data);

figure;
ka   = tiledlayout('flow'); title(ka  ,"Knee Flexion Angle");
figure;
kam  = tiledlayout('flow'); title(kam ,"Knee Adduction Moment");
figure;
kfm  = tiledlayout('flow'); title(kfm ,"Knee Flexion Moment");
figure;
ha   = tiledlayout('flow'); title(ha  ,"Hip Flexion Angle");
figure;
hfm  = tiledlayout('flow'); title(hfm ,"Hip Flexion Moment");
figure;
grfm = tiledlayout('flow'); title(grfm,"ML GRF");
figure;
grfa = tiledlayout('flow'); title(grfa,"AP GRF");
figure;
grfv = tiledlayout('flow'); title(grfv,"V GRF");
figure;
fa   = tiledlayout('flow'); title(fa  ,"Foot Segment Angle");

for c = 1:length(conditions)
    if isempty(data.(conditions{c})); continue; end
    nexttile(ka); hold on;
    make_plot(data.(conditions{c}).knee, conditions{c},      "angle",    1);
    nexttile(kam); hold on;
    make_plot(data.(conditions{c}).knee, conditions{c},     "moment",   2);
    nexttile(kfm); hold on;
    make_plot(data.(conditions{c}).knee, conditions{c},     "moment",   1);
    nexttile(ha); hold on;
    make_plot(data.(conditions{c}).hip, conditions{c},       "angle",    1);
    nexttile(hfm); hold on;
    make_plot(data.(conditions{c}).hip, conditions{c},      "moment",   1);
    nexttile(grfm); hold on;
    make_plot(data.(conditions{c}).grf, conditions{c},     "force",    1);
    nexttile(grfa); hold on;
    make_plot(data.(conditions{c}).grf, conditions{c},     "force",    2);
    nexttile(grfv); hold on;
    make_plot(data.(conditions{c}).grf, conditions{c},     "force",    3);
    nexttile(fa); hold on;
    make_plot(data.(conditions{c}).foot, conditions{c},      "angle",    1);

end


end

function make_plot(data2plot, condname,  type, axis)
    
rightraw_c = "#FF7174";
leftraw_c   = "#605F61";
right_c = "r";
left_c = "k";

subtitle(condname)
plot(squeeze(data2plot.r.(type).cycle_nd(:,axis,:)), Color=rightraw_c)
plot(squeeze(data2plot.l.(type).cycle_nd(:,axis,:)), Color=leftraw_c )
plot(data2plot.r.(type).cycle_nd_mean(:,axis), right_c, LineWidth=2, DisplayName="Right");
plot(data2plot.l.(type).cycle_nd_mean(:,axis), left_c , LineWidth=2, DisplayName="Left");
    
end
