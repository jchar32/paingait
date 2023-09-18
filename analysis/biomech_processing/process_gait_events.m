function [gait_events] = process_gait_events(all_data)
condition_names = fieldnames(all_data);
for c = 1:size(condition_names,1)
    
    % left 
    LON = all_data.(condition_names{c}).LON{:};
    LOFF = all_data.(condition_names{c}).LOFF{:};
    RON = all_data.(condition_names{c}).RON{:};
    ROFF = all_data.(condition_names{c}).ROFF{:};

    % gait events
    [gait_events.(condition_names{c})] = clean_gait_events(LON, LOFF, RON, ROFF);
end
end

function [events] = clean_gait_events(LON, LOFF, RON, ROFF)
    
    % remove any extra OFF events at start as needed (so first event is ON)
    LOFF = remove_extra_OFFs(LON, LOFF);
    ROFF = remove_extra_OFFs(RON, ROFF);
    
    % align events by removing errant ONs 
    [LON, LOFF] = align_events(LON, LOFF);
    [RON, ROFF] = align_events(RON, ROFF);

    assert(length(LON)==length(LOFF),"Event Timing needs manual intervention. There are mismatches between LONs and LOFFs")
    assert(length(RON)==length(ROFF),"Event Timing needs manual intervention. There are mismatches between RONs and ROFFs")
    
    events.r.ON = RON;
    events.l.ON = LON;
    events.r.OFF = ROFF;
    events.l.OFF = LOFF;
end

function OFF = remove_extra_OFFs(ON, OFF)
    isONfirst = ON(1) < OFF(1);
    while ~isONfirst
         OFF(1) = [];
         isONfirst = ON(1) < OFF(1);
    end
end

function [ON, OFF] = align_events(ON, OFF)
    % remove errant gait events so each ON and OFF event matches
    min_events = min([length(ON),length(OFF)]);
    timing_threshold = 1.8; 
    
    while length(ON) ~= length(OFF)
        ON(find(OFF(1:min_events)-ON(1:min_events) > timing_threshold, 1,"first")) = [];
    end
end


    
    