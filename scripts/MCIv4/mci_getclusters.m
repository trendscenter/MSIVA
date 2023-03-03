function D = mci_getclusters(B, minclustersize)

C = zeros(size(B)); C(find(B > 0)) = 1;
[L,NUM] = spm_bwlabel(C,26);
D = [];
cnum = 0;
for ii = 1:NUM
    IND = find(L == ii);
    numvox = length(IND);
    fprintf('Cluster %d, size = %d\n',ii,numvox )
    if numvox >= minclustersize;
        cnum = cnum+1;
        D.mask{cnum} = IND;
        D.n(cnum) = numvox;
    else

    end
end


