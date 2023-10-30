function pca_plot(X, Xmean, pc_percentiles_vec, axis)
    
    highlight_linewidth = 1.25;
    plot(X',Color="#d3d3d3", DisplayName="Raw") % raw waveforms
    plot(Xmean','k',LineWidth=2, DisplayName="Grand Mean") % mean across all waveforms
    
    % PC1 5th and 95th percentile
    plot(pc_percentiles_vec(1,:)', Color="#3960B4", LineWidth=highlight_linewidth, DisplayName="PC1 5th"); 
    plot(pc_percentiles_vec(2,:)', Color="#3960B4", LineWidth=highlight_linewidth,LineStyle="--", DisplayName="PC1 95th"); 
    % PC2 5th and 95th percentile
    if size(pc_percentiles_vec,1) > 2
    plot(pc_percentiles_vec(3,:)', Color="#B53AA3", LineWidth=highlight_linewidth, DisplayName="PC2 5th"); 
    plot(pc_percentiles_vec(4,:)', Color="#B53AA3", LineWidth=highlight_linewidth,LineStyle="--", DisplayName="PC2 95th");
    
    axes = ["Sagittal","Frontal","Transverse"];
    title([axes(axis) + " PCA"]);

    end

    if axis ==2
        f=get(gca,'Children');
        legend([f(6), f(5), f(4), f(3), f(2), f(1)], Location="northoutside",Orientation="horizontal")
    end
end