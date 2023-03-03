function updategradtype(O,gradt)

gradt = lower(gradt);
switch gradt
    case 'regular'
        O.gradtype = gradt;
    case 'relative'
        O.gradtype = gradt;
    case 'natural'
        O.gradtype = gradt;
        
    otherwise
        % put a warning here...
        O.gradtype = 'regular';
end

end
