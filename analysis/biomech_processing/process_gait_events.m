function [gait_events] = process_gait_events(all_data)
condition_names = fieldnames(all_data);
for c = 1:size(condition_names,1)
    if isempty(all_data.(condition_names{c})); continue; end
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
    [LON, LOFF, LON_next] = align_events(LON, LOFF);
    [RON, ROFF, RON_next] = align_events(RON, ROFF);

    assert(length(LON)==length(LOFF),"Event Timing needs manual intervention. There are mismatches between LONs and LOFFs")
    assert(length(RON)==length(ROFF),"Event Timing needs manual intervention. There are mismatches between RONs and ROFFs")
    
    events.r.ON = RON;
    events.l.ON = LON;
    events.r.OFF = ROFF;
    events.l.OFF = LOFF;
    events.r.ON_next = RON_next;
    events.l.ON_next = LON_next;
end

function OFF = remove_extra_OFFs(ON, OFF)
    isONfirst = ON(1) < OFF(1);
    while ~isONfirst
         OFF(1) = [];
         isONfirst = ON(1) < OFF(1);
    end
end

function [ON_cor, OFF_cor, ON_next_cor] = align_events(ON, OFF)
    % remove errant gait events so each ON and OFF event matches
    num_ONs = length(ON);
    num_OFFs = length(OFF);

    % is it a good stance?
    off_cor = nan(max(num_ONs, num_OFFs),1);
    on_cor = nan(max(num_ONs, num_OFFs),1);
    for i = 1:length(on_cor)
        temp = OFF - ON(i);
        possible_off=find(temp>0.2 & temp < 0.8);
        
        if isempty(possible_off)
        continue;
        elseif length(possible_off) > 2
           off_cor(i,1) = min(OFF(possible_off));
           on_cor(i,1) = ON(i); 
        else
        off_cor(i,1) = OFF(possible_off);
        on_cor(i,1) = ON(i);
        end
    end
    ON_cor = on_cor(~isnan(on_cor));
    OFF_cor = off_cor(~isnan(off_cor));
        
    % is it good cycle
    on_next_cor = nan(length(ON_cor),1);
    for i=1:length(ON_cor)-1
        temp2 = ON - ON_cor(i);
        possible_onnext = find(temp2 > 0.4);
        onnext_idx = ON(possible_onnext) > OFF_cor(i) & ON(possible_onnext) < OFF_cor(i+1);
        
        on_next_cor(i) = ON(possible_onnext(onnext_idx(1)));
    end
    ON_next_cor = on_next_cor(~isnan(on_next_cor));
end



    
    