function [waveforms] = compile_waveforms(all_data, events, sample_rate)

condition_names = fieldnames(all_data);
for c = 1:size(condition_names,1)
    
    waveforms.(condition_names{c}).hip.l.angle  = gather_strides(all_data.(condition_names{c}).LHA{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).hip.r.angle  = gather_strides(all_data.(condition_names{c}).RHA{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
    waveforms.(condition_names{c}).hip.l.moment = gather_strides(all_data.(condition_names{c}).LHM{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).hip.r.moment = gather_strides(all_data.(condition_names{c}).RHM{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
    
    waveforms.(condition_names{c}).knee.l.angle  = gather_strides(all_data.(condition_names{c}).LKA{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).knee.r.angle  = gather_strides(all_data.(condition_names{c}).RKA{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
    waveforms.(condition_names{c}).knee.l.moment = gather_strides(all_data.(condition_names{c}).LKM{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).knee.r.moment = gather_strides(all_data.(condition_names{c}).RKM{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
   
    waveforms.(condition_names{c}).ankle.l.angle  = gather_strides(all_data.(condition_names{c}).LAA{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).ankle.r.angle  = gather_strides(all_data.(condition_names{c}).RAA{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
    waveforms.(condition_names{c}).ankle.l.moment = gather_strides(all_data.(condition_names{c}).LAM{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).ankle.r.moment = gather_strides(all_data.(condition_names{c}).RAM{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);

    waveforms.(condition_names{c}).foot.l.angle  = gather_strides(all_data.(condition_names{c}).LF{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).foot.r.angle  = gather_strides(all_data.(condition_names{c}).RF{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
    
    waveforms.(condition_names{c}).shank.l.angle  = gather_strides(all_data.(condition_names{c}).LSK{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).shank.r.angle  = gather_strides(all_data.(condition_names{c}).RSK{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
    
    waveforms.(condition_names{c}).thigh.l.angle  = gather_strides(all_data.(condition_names{c}).LTH{1,1}, events.(condition_names{c}), "L", sample_rate.mocap);
    waveforms.(condition_names{c}).thigh.r.angle  = gather_strides(all_data.(condition_names{c}).RTH{1,1}, events.(condition_names{c}), "R", sample_rate.mocap);
  
    waveforms.(condition_names{c}).grf.l.force  = gather_strides(all_data.(condition_names{c}).FP1_filt{1,1}, events.(condition_names{c}), "L", sample_rate.analog);
    waveforms.(condition_names{c}).grf.r.force  = gather_strides(all_data.(condition_names{c}).FP2_filt{1,1}, events.(condition_names{c}), "R", sample_rate.analog);
  
end

end

function [out] = gather_strides(data, events, side, sample_rate)

    for s = 1:size(events.(lower(side)).ON)-1
        [ON, OFF, ON_next, stance_frames, stride_frames] = unpack_events(events, s, side, sample_rate);
        if ~(is_good_stance(ON*(1/sample_rate), OFF*(1/sample_rate)) && is_good_stride(ON*(1/sample_rate), ON_next*(1/sample_rate)))
            out.stance{s,1} = NaN;
            out.cycle{s,1} = NaN;
            out.stance_nd(:,1:3,s) = NaN;
            out.cycle_nd(:,1:3,s) = NaN;
            continue;
        end
        
        out.stance{s,1} = data(ON:OFF,:);
        out.cycle{s,1} = data(ON:ON_next,:);
        % normalize for 0-100%
        out.stance_nd(:,1:3,s) = time_normalize(data(ON:OFF,:));
        out.cycle_nd(:,1:3,s) = time_normalize(data(ON:ON_next,:));
    end

    out.stance_nd_mean = nanmean(out.stance_nd,3);
    out.stance_nd_sd = nanstd(out.stance_nd,[],3);

    out.cycle_nd_mean = nanmean(out.cycle_nd,3);
    out.cycle_nd_sd = nanstd(out.cycle_nd,[],3);

end

function [isgood]= is_good_stance(ON, OFF)
    isgood=true;
    if (OFF - ON) > 1.5
        isgood=false;
    elseif (OFF - ON) < 0
        isgood=false;
    end
end

function [isgood] = is_good_stride(ON, ON_next)
    isgood=true;
    if (ON_next - ON) > 2
        isgood=false;
    elseif (ON_next - ON) < 0
         isgood=false;
    end
end

function [nd] = time_normalize(data)
    nd=interp1(1:size(data,1), data, linspace(0, size(data,1), 100 ),"pchip");
end

