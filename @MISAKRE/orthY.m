function Vt = orthY(O)

Vt = cell(size(O.W));
for mm = O.M
    [t1, D] = eig(O.Y{mm}*O.Y{mm}');
    D = diag(D);
    Vt{mm} = bsxfun(@times,1./D,t1')*O.Y{mm};
end

end