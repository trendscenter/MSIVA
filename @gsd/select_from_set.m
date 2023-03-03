function [selected, done] = select_from_set(O, set_, how_many, skip)
    selected = set_(1:skip:end);
    if length(selected) < how_many
        done = false;
    else
        selected = selected(1:how_many);
        done = true;
    end
end