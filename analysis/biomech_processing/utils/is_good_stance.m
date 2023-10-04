function [isgood]= is_good_stance(ON, OFF)
    isgood=true;
    if (OFF - ON) > 0.8
        isgood=false;
    elseif (OFF - ON) < 0
        isgood=false;
    end
end

