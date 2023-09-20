function [discrete_data] = calculate_discrete_outcomes(all_data, events, sample_rate)

condition_names = fieldnames(all_data);
for c = 1:1%size(condition_names,1)
    
    % left 
    LON =   events.(condition_names{c}).r.ON(:);
    LOFF =  events.(condition_names{c}).l.OFF(:);
    RON =   events.(condition_names{c}).r.ON(:);
    ROFF =  events.(condition_names{c}).l.OFF(:);

     % spatiotemporal
    [discrete_data.(condition_names{c}).r.temporal] = temporal_outcomes(RON, ROFF);
    [discrete_data.(condition_names{c}).l.temporal] = temporal_outcomes(LON, LOFF);

    % Kinematic and Kinetic Outcomes    
    for s = 1:size(events.(condition_names{c}).l.ON)-1
        discrete_data.(condition_names{c}).knee.l = knee_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), s, "L", sample_rate);
        discrete_data.(condition_names{c}).ankle.l = ankle_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), s, "L", sample_rate);
    end
    
    for s = 1:size(events.(condition_names{c}).r.ON)-1
        discrete_data.(condition_names{c}).knee.r = knee_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), s, "R", sample_rate);
        discrete_data.(condition_names{c}).ankle.r = ankle_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), s, "R", sample_rate);
    end
    

end


end


function [out] = knee_outcomes(data, events, stride, side, sample_rate)

[ON, OFF, ON_next, stance_frames, stride_frames] = unpack_events(events, stride, side, sample_rate);

% Kinematic  outcomes
KA = data.(side + "KA"){1,1};
% Sagittal
[out.peak_fa,      out.time.peak_fa]           = max(KA(ON:ON+round(stance_frames*0.5), 1));
[out.peak_ea,      out.time.peak_ea]           = min(KA(out.time.peak_fa:OFF, 1)); out.time.peak_ea + [0, + out.time.peak_fa]; 
[out.peak_fa_swing,out.time.peak_fa_swing]     = max(KA(ON:ON_next, 1));
[out.fa_hs,    out.time.fa_hs]         = max(KA(ON, 1));
[out.fa_to,    out.time.fa_to]         = max(KA(OFF, 1));

[out.peak_fa_swing,out.time.peak_fa_swing]     = max(KA(ON:ON_next, 1));
[out.fa_hs,    out.time.fa_hs]         = max(KA(ON, 1));
[out.fa_to,    out.time.fa_to]         = max(KA(OFF, 1));

% Frontal
[out.peak_aa,      out.time.peak_aa]           = max(KA(ON:OFF, 2));
[out.mean_aa]           = mean(KA(ON:OFF, 2));

[out.peak_aa_swing,out.time.peak_aa_swing]     = max(KA(ON:ON_next, 2));
[out.aa_hs,    out.time.aa_hs]         = max(KA(ON, 2));
[out.aa_to,    out.time.aa_to]         = max(KA(OFF, 2));

% angular velocity
[ka_deriv] = central_difference(KA(:,:), sample_rate.mocap);


% Kinetic outcomes
KM = data.(side + "KM"){1,1};
% Sagittal
[out.peak_fm, out.time.peak_fm] = max(KM(ON:OFF, 1));
[out.peak_em, out.time.peak_em] = min(KM(ON:OFF, 1));
positive_idx = find(KM(ON:OFF,1) >= 0) + ON;
[out.impulse_fm]                    = trapz(KM(positive_idx,1))*(1/sample_rate.mocap);

% Frontal
[out.peak_am, out.time.peak_am]     = max(KM(ON:OFF,2));
[out.peak_am1, out.time.peak_am1]   = max(KM(ON:ON+round(stance_frames*0.5), 2));
[out.peak_am2, out.time.peak_am2]   = max(KM(ON+round(stance_frames*0.5):OFF, 2)); out.time.peak_am2 + ON+round(stance_frames*0.5);
[out.peak_am, out.time.peak_am]     = min(KM(out.time.peak_am1:out.time.peak_am2, 2)); out.time.peak_am + out.time.peak_am1;
positive_idx = find(KM(ON:OFF,2) >= 0) + ON;
[out.impulse_am]                    = trapz(KM(positive_idx,2))*(1/sample_rate.mocap);

% loading rate
[km_deriv] = central_difference(KM(:,:), sample_rate.mocap);
[out.peak_fm_rate, out.time.peak_fm_rate] = max(km_deriv(ON:ON+round(stance_frames*0.5),1));
[out.peak_am_rate, out.time.peak_am_rate] = min(km_deriv(ON:ON+round(stance_frames*0.5),2));

% Composite outcomes
% per Zeni et al 2009 Clin Biomech (24) - 3% stance to to peak knee flexion angle
[out.djs_f, out.djs_a] = dynamic_joint_stiffness(km_deriv, ka_deriv, round(stance_frames*0.03) + ON, out.time.peak_fa+ON);


end



function [out] = ankle_outcomes(data, events, stride, side, sample_rate)
[ON, OFF, ON_next, stance_frames, stride_frames] = unpack_events(events, stride, side, sample_rate);

% Kinematic Outcomes
AA = data.(side + "AA"){1,1};
% Sagittal
[out.peak_dfa,      out.time.peak_dfa]           = min(AA(ON:ON+round(stance_frames*0.5), 1));
[out.peak_pfa,      out.time.peak_pfa]           = min(AA(ON+round(stance_frames*0.5):OFF, 1));
[out.fa_hs,      out.time.fa_hs]           = max(AA(ON,1));
[out.fa_to,      out.time.fa_to]           = max(AA(OFF,1));

% Frontal
[out.peak_eva,      out.time.peak_eva]           = min(AA(ON:OFF, 2));
[out.eva_hs,      out.time.eva_hs]           = max(AA(ON,2));
[out.eva_to,      out.time.eva_to]           = max(AA(OFF,2));

% Angular velocity
[aa_deriv] = central_difference(AA(:,:), sample_rate.mocap);

% Kinetic
AM = data.(side + "AM"){1,1};
[out.peak_dfm,      out.time.peak_dfm]           = max(AM(ON:OFF, 1));
[out.peak_pfm,      out.time.peak_pfm]           = min(AM(ON:ON+round(stance_frames*0.25), 1));
[out.peak_evm,      out.time.peak_evm]           = min(AM(ON:OFF, 2));

% Loading rate
[am_deriv] = central_difference(AM(:,:), sample_rate.mocap);

% Composite
[out.djs_pf, out.djs_a] = dynamic_joint_stiffness(am_deriv, aa_deriv, round(stance_frames*0.03) + ON, out.time.peak_eva+ON);

end



function [events_out] = temporal_outcomes(ON, OFF)
    
    events_out.stride_time = diff(ON);
    events_out.stride_time(events_out.stride_time > 1.8) = NaN; % replace the swing times that are incorrect due to missing ONs that were removed
    events_out.stance_time = OFF-ON; 
    events_out.swing_time = ON(2:end) - OFF(1:end-1); 
    events_out.swing_time(events_out.swing_time > 1.2) = NaN; % replace the swing times that are incorrect due to missing ONs that were removed
    events_out.stance_time_perc = events_out.stance_time(1:end-1) ./ events_out.stride_time;
    events_out.swing_time_perc = events_out.swing_time(1:end) ./ events_out.stride_time;
    events_out.cadence = sum(~isnan(events_out.stride_time))/(nansum(events_out.stride_time)/60);
end

function [isgood]= is_good_stride(ON, OFF,s)
    isgood=true;
    if (OFF(s) - ON(s)) > 1.5
        isgood=false;
    elseif (OFF(s) - ON(s) < 0)
    end

end

function [dy] = central_difference(y, h)
% expects y is an n samples by m channels array
n =size(y,1);
for i=2:n-1
    dy(i,:) = (y(i-1,:)+y(i+1,:)) ./2 ./h;
end
end

function [ON, OFF, ON_next, stance_frames, stride_frames] = unpack_events(events, stride, side, sample_rate)
% unpack events
ON = round(events.(lower(side)).ON(stride) * sample_rate.mocap);
OFF = round(events.(lower(side)).OFF(stride) * sample_rate.mocap);
ON_next = round(events.(lower(side)).ON(stride+1) * sample_rate.mocap);
stance_frames = OFF-ON;
stride_frames = ON_next-ON;
end

function [djs_f, djs_a] = dynamic_joint_stiffness(moment_deriv, angle_deriv, startframe, endframe)
% Composite outcomes
% per Zeni et al 2009 Clin Biomech (24) - 3% stance to to peak knee flexion angle

[fp,~]=polyfit(angle_deriv(startframe:endframe,1), moment_deriv(startframe:endframe,1), 1);
[ap,~]=polyfit(angle_deriv(startframe:endframe,2), moment_deriv(startframe:endframe,1), 2);
[djs_f] = fp(1); % slope of line
[djs_a] = ap(1); 


end