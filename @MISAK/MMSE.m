function [mMSE] = MMSE(O,Y)

mMSE = O.ut.MMSE(Y,O.Y,O.S);

end