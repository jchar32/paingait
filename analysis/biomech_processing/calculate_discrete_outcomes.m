function [discrete_data] = calculate_discrete_outcomes(all_data, events, sample_rate)

condition_names = fieldnames(all_data);
for c = 1:size(condition_names,1)
    if isempty(all_data.(condition_names{c})); continue; end
    % Events 
    LON =   events.(condition_names{c}).l.ON(:);
    LOFF =  events.(condition_names{c}).l.OFF(:);
    LON_next = events.(condition_names{c}).l.ON_next(:);
    RON =   events.(condition_names{c}).r.ON(:);
    ROFF =  events.(condition_names{c}).r.OFF(:);
    RON_next = events.(condition_names{c}).r.ON_next(:);

    % spatiotemporal
    [discrete_data.(condition_names{c}).temporal.r] = temporal_outcomes(RON, ROFF);
    [discrete_data.(condition_names{c}).temporal.l] = temporal_outcomes(LON, LOFF);
    
    % counter to track number of good strides that are used in discrete outcomes
    discrete_data.(condition_names{c}).left_stride_counter = 0;
    discrete_data.(condition_names{c}).right_stride_counter = 0;
    
    % Kinematic and Kinetic Outcomes    
%     for s = 1:size(events.(condition_names{c}).l.ON)-1 % Left Leg

%         [ON, OFF, ON_next, ~, ~] = unpack_events(events.(condition_names{c}), s, "L", sample_rate.mocap);
%         if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
%             continue;
%         end
        discrete_data.(condition_names{c}).left_stride_counter =  discrete_data.(condition_names{c}).left_stride_counter+1;

        discrete_data.(condition_names{c}).knee.l = knee_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "L", sample_rate);
        discrete_data.(condition_names{c}).ankle.l = ankle_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "L", sample_rate);
        discrete_data.(condition_names{c}).hip.l = hip_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), "L", sample_rate);
        
        discrete_data.(condition_names{c}).grf.l = grf_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "L", sample_rate);

        discrete_data.(condition_names{c}).foot.l = foot_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "L", sample_rate);
        discrete_data.(condition_names{c}).shank.l = shank_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "L", sample_rate);
        discrete_data.(condition_names{c}).thigh.l = thigh_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "L", sample_rate);
%     end
    
%     for s = 1:size(events.(condition_names{c}).r.ON)-1 % Right Leg

%         [ON, OFF, ON_next, ~, ~] = unpack_events(events.(condition_names{c}), s, "R", sample_rate.mocap);
%         if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
%             continue;
%         end
        discrete_data.(condition_names{c}).right_stride_counter =  discrete_data.(condition_names{c}).right_stride_counter +1;

        discrete_data.(condition_names{c}).knee.r = knee_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), "R", sample_rate);
        discrete_data.(condition_names{c}).ankle.r = ankle_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "R", sample_rate);
        discrete_data.(condition_names{c}).hip.r = hip_outcomes(all_data.(condition_names{c}), events.(condition_names{c}), "R", sample_rate);
        
        discrete_data.(condition_names{c}).grf.r = grf_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "R", sample_rate);
        
        discrete_data.(condition_names{c}).foot.r = foot_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "R", sample_rate);
        discrete_data.(condition_names{c}).shank.r = shank_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "R", sample_rate);
        discrete_data.(condition_names{c}).thigh.r = thigh_outcomes(all_data.(condition_names{c}), events.(condition_names{c}),  "R", sample_rate);

%     end
end
end

function [out] = hip_outcomes(data, events, side, sample_rate)

    
    for stride = 1:size(events.(lower(side)).ON)-1 % Right Leg
    
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
    
        HA = data.(side + "HA"){1,1};
        % Sagittal
        [out.peak_fa(stride),      out.time.peak_fa(stride)]           = max(HA(ON:ON+round(stance_frames*0.5), 1));
        [out.peak_ea(stride),      out.time.peak_ea(stride)]           = min(HA(ON:OFF, 1));
        [out.peak_fa_swing(stride), out.time.peak_fa_swing(stride)]     = max(HA(ON:ON_next, 1));
        [out.fa_hs(stride),    out.time.fa_hs(stride)]         = max(HA(ON, 1));
        [out.fa_to(stride),    out.time.fa_to{stride}]         = max(HA(OFF, 1));
    
        % Frontal
        [out.peak_aa(stride),       out.time.peak_aa(stride)]           = min(HA(ON:OFF, 2));
        [out.mean_aa(stride)]           = mean(HA(ON:OFF, 2));
        [out.peak_aa_swing(stride), out.time.peak_aa_swing(stride)]     = max(HA(ON:ON_next, 2));
        [out.aa_hs(stride),         out.time.aa_hs(stride)]         = max(HA(ON, 2));
        [out.aa_to(stride),         out.time.aa_to(stride)]         = max(HA(OFF, 2));
    
        % angular velocity
        [ha_deriv] = central_difference(HA(:,:), sample_rate.mocap);
    
        % Kinetics
        HM = data.(side + "HM"){1,1};
        % Sagittal
        [out.peak_fm(stride), out.time.peak_fm(stride)] = min(HM(ON:OFF, 1));
        [out.peak_em(stride), out.time.peak_em(stride)] = max(HM(ON:OFF, 1));
        positive_idx = find(HM(ON:OFF,1) >= 0) + ON;
        neg_idx = find(HM(ON:OFF,1) <= 0) + ON;
        [out.impulse_em(stride)]                    = trapz(HM(positive_idx,1))*(1/sample_rate.mocap);
        [out.impulse_fm(stride)]                    = trapz(HM(neg_idx,1))*(1/sample_rate.mocap);
    
        % Frontal
        [out.peak_am(stride), out.time.peak_am(stride)]     = max(HM(ON:OFF,2));
        [out.peak_am1(stride), out.time.peak_am1(stride)]   = max(HM(ON:ON+round(stance_frames*0.5), 2));
        [out.peak_am2(stride), out.time.peak_am2(stride)]   = max(HM(ON+round(stance_frames*0.5):OFF, 2));
        out.time.peak_am2(stride) = out.time.peak_am2(stride) + ON+round(stance_frames*0.5);
        [out.unload_am(stride), out.time.unload_am(stride)]     = min(HM(out.time.peak_am1(stride):out.time.peak_am2(stride), 2));
        out.time.unload_am(stride) = out.time.unload_am(stride) + out.time.peak_am1(stride);
        positive_idx = find(HM(ON:OFF,2) >= 0) + ON;
        [out.impulse_am(stride)]                    = trapz(HM(positive_idx,2))*(1/sample_rate.mocap);
    
        % loading rate
        [hm_deriv] = central_difference(HM(:,:), sample_rate.mocap);
        [out.peak_fm_rate(stride), out.time.peak_fm_rate(stride)] = max(hm_deriv(ON:ON+round(stance_frames*0.5),1));
        [out.peak_am_rate(stride), out.time.peak_am_rate(stride)] = min(hm_deriv(ON:ON+round(stance_frames*0.5),2));
    
        % Composite outcomes
        % per Zeni et al 2009 Clin Biomech (24) - 3% stance to to peak knee flexion angle
        % [out.djs_f, out.djs_a] = dynamic_joint_stiffness(hm_deriv, ha_deriv, round(stance_frames*0.03) + ON, out.time.peak_fa+ON);
    end
end

function [out] = knee_outcomes(data, events, side, sample_rate)


    for stride = 1:size(events.(lower(side)).ON)-1 % Right Leg
    
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
    
        % Kinematic  outcomes
        KA = data.(side + "KA"){1,1};
        % Sagittal
        [out.peak_fa(stride),      out.time.peak_fa(stride)]           = max(KA(ON:ON+round(stance_frames*0.5), 1));
        [out.peak_ea(stride),      out.time.peak_ea(stride)]           = min(KA(out.time.peak_fa:OFF, 1));
        out.time.peak_ea(stride) = out.time.peak_ea(stride) + out.time.peak_fa(stride);
        [out.peak_fa_swing(stride), out.time.peak_fa_swing(stride)]     = max(KA(ON:ON_next, 1));
        [out.fa_hs(stride),    out.time.fa_hs(stride)]         = max(KA(ON, 1));
        [out.fa_to(stride),    out.time.fa_to(stride)]         = max(KA(OFF, 1));
    
        % Frontal
        [out.peak_aa(stride),      out.time.peak_aa(stride)]           = max(KA(ON:OFF, 2));
        [out.mean_aa(stride)]           = mean(KA(ON:OFF, 2));
    
        [out.peak_aa_swing(stride), out.time.peak_aa_swing(stride)]     = max(KA(ON:ON_next, 2));
        [out.aa_hs(stride),    out.time.aa_hs(stride)]         = max(KA(ON, 2));
        [out.aa_to(stride),    out.time.aa_to(stride)]         = max(KA(OFF, 2));
    
        % angular velocity
        [ka_deriv] = central_difference(KA(:,:), sample_rate.mocap);
    
    
        % Kinetic outcomes
        KM = data.(side + "KM"){1,1};
        % Sagittal
        [out.peak_fm(stride), out.time.peak_fm(stride)] = max(KM(ON:OFF, 1));
        [out.peak_em(stride), out.time.peak_em(stride)] = min(KM(ON:OFF, 1));
        positive_idx = find(KM(ON:OFF,1) >= 0) + ON;
        [out.impulse_fm(stride)]                    = trapz(KM(positive_idx,1))*(1/sample_rate.mocap);
    
        % Frontal
        [out.peak_am(stride), out.time.peak_am(stride)]     = max(KM(ON:OFF,2));
        [out.peak_am1(stride), out.time.peak_am1(stride)]   = max(KM(ON:ON+round(stance_frames*0.5), 2));
        [out.peak_am2(stride), out.time.peak_am2(stride)]   = max(KM(ON+round(stance_frames*0.5):OFF, 2));
        out.time.peak_am2(stride)= out.time.peak_am2(stride) + ON+round(stance_frames*0.5);
        [out.unload_am(stride), out.time.unload_am(stride)]     = min(KM(out.time.peak_am1(stride):out.time.peak_am2(stride), 2));
        out.time.unload_am(stride)= out.time.unload_am(stride) + out.time.peak_am1(stride);
        positive_idx = find(KM(ON:OFF,2) >= 0) + ON;
        [out.impulse_am(stride)]                    = trapz(KM(positive_idx,2))*(1/sample_rate.mocap);
    
        % loading rate
        [km_deriv] = central_difference(KM(:,:), sample_rate.mocap);
        [out.peak_fm_rate(stride), out.time.peak_fm_rate(stride)] = max(km_deriv(ON:ON+round(stance_frames*0.5),1));
        [out.peak_am_rate(stride), out.time.peak_am_rate(stride)] = min(km_deriv(ON:ON+round(stance_frames*0.5),2));
    
        % Composite outcomes
        % per Zeni et al 2009 Clin Biomech (24) - 3% stance to to peak knee flexion angle
        [out.djs_f(stride), out.djs_a(stride)] = dynamic_joint_stiffness(km_deriv, ka_deriv, round(stance_frames*0.03) + ON, out.time.peak_fa(stride)+ON);
    end
end

function [out] = ankle_outcomes(data, events, side, sample_rate)

    for stride = 1:size(events.(lower(side)).ON)-1
    
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
    
        % Kinematic Outcomes
        AA = data.(side + "AA"){1,1};
        % Sagittal
        [out.peak_dfa(stride),      out.time.peak_dfa(stride)]           = min(AA(ON:ON+round(stance_frames*0.5), 1));
        [out.peak_pfa(stride),      out.time.peak_pfa(stride)]           = min(AA(ON+round(stance_frames*0.5):OFF, 1));
        [out.fa_hs(stride),      out.time.fa_hs(stride)]           = max(AA(ON,1));
        [out.fa_to(stride),      out.time.fa_to(stride)]           = max(AA(OFF,1));
    
        % Frontal
        [out.peak_eva(stride),      out.time.peak_eva(stride)]           = min(AA(ON:OFF, 2));
        [out.eva_hs(stride),      out.time.eva_hs(stride)]           = max(AA(ON,2));
        [out.eva_to(stride),      out.time.eva_to(stride)]           = max(AA(OFF,2));
    
        % Angular velocity
        [aa_deriv] = central_difference(AA(:,:), sample_rate.mocap);
    
        % Kinetic
        AM = data.(side + "AM"){1,1};
        [out.peak_dfm(stride),      out.time.peak_dfm(stride)]           = max(AM(ON:OFF, 1));
        [out.peak_pfm(stride),      out.time.peak_pfm(stride)]           = min(AM(ON:ON+round(stance_frames*0.25), 1));
        [out.peak_evm(stride),      out.time.peak_evm(stride)]           = min(AM(ON:OFF, 2));
    
        % Loading rate
        [am_deriv] = central_difference(AM(:,:), sample_rate.mocap);
    
        % Composite
        % [out.djs_pf, out.djs_a] = dynamic_joint_stiffness(am_deriv, aa_deriv, round(stance_frames*0.03) + ON, out.time.peak_eva+ON);
    end
end

function [out] = grf_outcomes(data, events, side, sample_rate)

    for stride = 1:size(events.(lower(side)).ON)-1 % Right Leg
    
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
    
        % GRF signal order = [ML,AP,V];
        % GRF convention: + = medial, propulsive/forward, up
    
        if side == "L"
            GRF = data.FP1_filt{1,1};
        else
            GRF = data.FP2_filt{1,1} .* [-1,1,1];
        end
    
        % Mediolateral
        [out.peak_ml(stride), out.time.peak_ml(stride)]     = max(GRF(ON:OFF,1));
        [out.peak_ml1(stride), out.time.peak_ml1(stride)]   = max(GRF(ON:ON+round(stance_frames*0.5), 1));
        [out.peak_ml2(stride), out.time.peak_ml2(stride)]   = max(GRF(ON+round(stance_frames*0.5):OFF, 1));
        out.time.peak_ml2(stride) = out.time.peak_ml2(stride) + ON+round(stance_frames*0.5);
        [out.unload_ml(stride), out.time.unload_ml(stride)]     = min(GRF(out.time.peak_ml1:out.time.peak_ml2, 1));
        out.time.unload_ml(stride)= out.time.unload_ml(stride) + out.time.peak_ml1(stride);
        positive_idx = find(GRF(ON:OFF,1) >= 0) + ON;
        [out.impulse_ml(stride)]                    = trapz(GRF(positive_idx,1))*(1/sample_rate.analog);
    
        % Anteroposterior
        [out.peak_brake(stride), out.time.peak_brake(stride)] = min(GRF(ON:ON+round(stance_frames*0.5), 2));
        [out.peak_prop(stride), out.time.peak_prop(stride)] = max(GRF(ON:OFF, 2));
        neg_idx = find(GRF(ON:OFF,2) <= 0) + ON;
        [out.impulse_brake(stride)]                    = trapz(GRF(neg_idx,2))*(1/sample_rate.analog);
        positive_idx = find(GRF(ON:OFF,2) >= 0) + ON;
        [out.impulse_prop(stride)]                    = trapz(GRF(positive_idx,2))*(1/sample_rate.analog);
    
        % Vertical
        [out.peak_v(stride), out.time.peak_v(stride)]     = max(GRF(ON:OFF,3));
        [out.peak_v1(stride), out.time.peak_v1(stride)]   = max(GRF(ON:ON+round(stance_frames*0.5), 3));
        [out.peak_v2(stride), out.time.peak_v2(stride)]   = max(GRF(ON+round(stance_frames*0.5):OFF, 3));
        out.time.peak_v2(stride) = out.time.peak_v2(stride) + ON+round(stance_frames*0.5);
        [out.unload_v(stride), out.time.unload_v(stride)]     = min(GRF(out.time.peak_v1(stride):out.time.peak_v2(stride), 3));
        out.time.unload_v(stride)= out.time.unload_v(stride) + out.time.peak_v1(stride);
        positive_idx = find(GRF(ON:OFF,3) >= 0) + ON;
        [out.impulse_v(stride)]                    = trapz(GRF(positive_idx,3))*(1/sample_rate.analog);
    
        % loading rate
        [grf_deriv] = central_difference(GRF(:,:), sample_rate.analog);
        [out.peak_ml_rate(stride), out.time.peak_ml_rate(stride)] = max(grf_deriv(ON:ON+round(stance_frames*0.5),1));
        [out.peak_ap_rate(stride), out.time.peak_ap_rate(stride)] = min(grf_deriv(ON:ON+round(stance_frames*0.5),2));
        [out.peak_v_rate(stride), out.time.peak_v_rate(stride)] = min(grf_deriv(ON:ON+round(stance_frames*0.5),3));
    
    end
end

function [out] = foot_outcomes(data, events, side, sample_rate)

    for stride = 1:size(events.(lower(side)).ON)-1 % Right Leg
    
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
        
        % Kinematic Outcomes
        FA = data.(side + "F"){1,1};
        % not actually flexion(fa)/extension angles(ea), but this is fine for now.
        [out.peak_pfa(stride),      out.time.peak_pfa(stride)]           = max(FA(ON:OFF,1));
        [out.peak_dfa(stride),      out.time.peak_dfa(stride)]           = min(FA(ON:OFF,1));
        
        [out.a_hs(stride),      out.time.a_hs(stride)]           = max(FA(ON,1));
        [out.a_to(stride),      out.time.a_to(stride)]           = max(FA(OFF,1));
        [out.a_ms(stride),      out.time.a_ms(stride)]           = max(FA(ON+round(stance_frames*0.5),1)) ;
        out.time.a_ms(stride) = out.time.a_ms(stride) + ON;
        [out.a_termswing(stride)]                        = max(FA(ON-round(stance_frames*0.1):ON,1));
    end
end

function [out] = shank_outcomes(data, events, side, sample_rate)

    for stride = 1:size(events.(lower(side)).ON)-1 % Right Leg
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
        
        % Kinematic Outcomes
        SA = data.(side + "SK"){1,1};
        if  side == "R"
            SA = SA.*[1,1,-1];
        end
        % not actually flexion(fa)/extension angles(ea), but this is fine for now.
        [out.peak_fa(stride),      out.time.peak_fa(stride)]           = max(SA(ON:OFF,1));
        [out.peak_ea(stride),      out.time.peak_ea(stride)]           = min(SA(ON:OFF,1));
        
        [out.a_hs(stride),      out.time.a_hs(stride)]           = max(SA(ON,1));
        [out.a_to(stride),      out.time.a_to(stride)]           = max(SA(OFF,1));
        [out.a_ms(stride),      out.time.a_ms(stride)]           = max(SA(ON+round(stance_frames*0.5),1)) ;
        out.time.a_ms(stride) = out.time.a_ms(stride) + ON;
        [out.a_termswing(stride)]                        = max(SA(ON-round(stance_frames*0.1):ON,1));
    end
end

function [out] = thigh_outcomes(data, events, side, sample_rate)
    
    for stride = 1:size(events.(lower(side)).ON)-1 % Right Leg
    
        [ON, OFF, ON_next, stance_frames, ~] = unpack_events(events, stride, side, sample_rate.mocap);
        if ~(is_good_stance(ON*(1/sample_rate.mocap), OFF*(1/sample_rate.mocap)) && is_good_stride(ON*(1/sample_rate.mocap), ON_next*(1/sample_rate.mocap)))
            continue;
        end
    
        % Kinematic Outcomes
        TA = data.(side + "TH"){1,1};
        if  side == "R"
            TA = TA.*[1,1,-1];
        end
        % not actually flexion(fa)/extension angles(ea), but this is fine for now.
        [out.peak_fa(stride),      out.time.peak_fa(stride)]           = max(TA(ON:OFF,1));
        [out.peak_ea(stride),      out.time.peak_ea(stride)]           = min(TA(ON:OFF,1));
    
        [out.a_hs(stride),      out.time.a_hs(stride)]           = max(TA(ON,1));
        [out.a_to(stride),      out.time.a_to(stride)]           = max(TA(OFF,1));
        [out.a_ms(stride),      out.time.a_ms(stride)]           = max(TA(ON+round(stance_frames*0.5),1)) ;
        out.time.a_ms(stride) = out.time.a_ms(stride) + ON;
        [out.a_termswing(stride)]                        = max(TA(ON-round(stance_frames*0.1):ON,1));
    end
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


%% Helper functions
function [dy] = central_difference(y, fs)
    % y = n x m where n=samples and m=channels
    % fs = sampling freq
    n = size(y,1);
    dy = nan(size(y));
    dt=1/fs;

    % forward diff
    dy(1,:) = (y(2,:) - y(1,:)) / dt;

    % backward diff
    dy(n,:) = (y(n,:) - y(n-1,:)) / dt;
    
    % central
    dy(2:n-1,:) = (y(3:n,:) - y(1:n-2,:)) ./ (2*dt);

end

function [djs_f, djs_a] = dynamic_joint_stiffness(moment_deriv, angle_deriv, startframe, endframe)
    % Composite outcomes
    % per Zeni et al 2009 Clin Biomech (24) - 3% stance to to peak knee flexion angle
    
    [fp,~]=polyfit(angle_deriv(startframe:endframe,1), moment_deriv(startframe:endframe,1), 1);
    [ap,~]=polyfit(angle_deriv(startframe:endframe,2), moment_deriv(startframe:endframe,1), 2);
    [djs_f] = fp(1); 
    [djs_a] = ap(1); 
end