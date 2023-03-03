function setREapproach(O,REa)

switch upper(REa)
    case 'PINV'
        O.REapproach_ = true;
        O.REapproach = 'PINV';
    case 'WT'
        O.REapproach_ = false;
        O.REapproach = 'WT';
        
    otherwise
        % put a warning here...
        O.REapproach_ = true;
        O.REapproach = 'PINV';
end

end