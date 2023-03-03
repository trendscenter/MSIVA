function saveX(O,path,prefix,suffix)

X = O.genX();
save(fullfile(path,[prefix '_' O.tag '_data_' suffix '.mat']), 'X');

end