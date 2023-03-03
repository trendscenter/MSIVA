function W = unstackW(w,M,C,V)

% Unstack wout: *wout is a col. vector. *Number of elements in Wf = Cf * Nf, different for each f
numel_Wf = C(M).*V(M);
% To partition wout into M Wf's, need start/end indexes corresponding to each f.
ixe = cumsum(numel_Wf);       % End index
ixb = 1 + [0 ixe(1:(end-1))]; % Begin index
% Parse wout:  Elements of wout are rows of Wf (M=1 first).
W = cell(1,length(M));
for mm = M
    W{mm} = reshape( w(ixb(M==mm):ixe(M==mm)), V(mm), C(mm) )';
end
W = W(M);