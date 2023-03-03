function setREtype(O,REt)

switch upper(REt)
    case 'NMSE'
        O.REtype_ = true;
        O.REtype = 'NMSE';
    case 'MSE'
        O.REtype_ = false;
        O.REtype = 'MSE';
    otherwise
        O.REtype_ = {};
end

end