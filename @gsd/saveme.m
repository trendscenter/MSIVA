function saveme(O,path,prefix,suffix)

save(fullfile(path,[prefix '_' O.tag '_GT_' suffix '.mat']), 'O');

end