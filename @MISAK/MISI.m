function [mISI] = MISI(O,varargin)

if ~isempty(varargin)
    var = varargin{1,1};
    A = var{1,1};
    Sgt = var{1,2};
end

S = cellfun(@(s) full(s), Sgt, 'Un', 0);
mISI = O.ut.MISI(O.W,A,S);

end