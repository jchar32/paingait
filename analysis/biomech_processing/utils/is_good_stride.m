function [isgood] = is_good_stride(ON, ON_next)
    isgood=true;
    if (ON_next - ON) > 1.5
        isgood=false;
    elseif (ON_next - ON) < 0
        isgood=false;
    end
end
