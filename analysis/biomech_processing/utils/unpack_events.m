function [ON, OFF, ON_next, stance_frames, stride_frames] = unpack_events(events, stride, side, sample_rate)
    % unpack events
    ON = round(events.(lower(side)).ON(stride) * sample_rate);
    OFF = round(events.(lower(side)).OFF(stride) * sample_rate);
    ON_next = round(events.(lower(side)).ON(stride+1) * sample_rate);
    stance_frames = OFF-ON;
    stride_frames = ON_next-ON;
end

