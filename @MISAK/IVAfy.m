function Sold = IVAfy(O)

Sold = O.S;

S_ = cell(size(O.S));
for mm = O.M
    ixs = [];
    for kk = find(O.nes)'
        ixs = [ixs find(O.S{mm}(kk,:) == 1, 1)];
    end
    S_{mm} = O.S{mm}(:,ixs);
end
O.updateCS(S_); % Update S and C

end
