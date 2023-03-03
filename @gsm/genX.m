function X = genX(O)

X = cellfun(@(a,y,n) a*y + n, O.A, O.Y, O.noise, 'Un', 0);

end