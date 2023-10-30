function mean_wfrm_plot(Xmean, condmeanstbl, axis)

% 10 color palette repeated twice
colors=["#f94144", "#f3722c", "#f8961e", "#f9844a", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1", ...
    "#f94144", "#f3722c", "#f8961e", "#f9844a", "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1"]; 
    
axes = ["Sagittal","Frontal","Transverse"];
    varnames = fieldnames(condmeanstbl);
    colors2use = colors(1:1:length(varnames)-3);
    for i=1:length(varnames)-3
        plot(condmeanstbl.(varnames{i}), LineWidth=1.25, Color=colors2use(i));
    end
    plot(Xmean', 'k', LineWidth=1.25)
    
    if axis ==2
        legend(varnames(1:end-3), Location="northoutside", Orientation="horizontal")
    end
    title([axes(axis) + " Mean Waveforms"])
    xlabel("Percent Gait Cycle (%)")
  
end