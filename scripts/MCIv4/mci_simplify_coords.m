function C = mci_simplify_coords(coords)

C.x = squeeze(coords(1,:,1,1));
C.y = squeeze(coords(2,1,:,1))';
C.z = squeeze(coords(3,1,1,:))';