function [mISI] = MISI(O,A)

S = cellfun(@(s) full(s), O.S, 'Un', 0);
mISI = O.ut.MISI(O.W,A,S);

end