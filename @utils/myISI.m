function ISI = myISI(WA)

N = length(WA);

WA = abs(WA);

ISI = 0;
ISI = ISI + sum(sum(WA,2)./max(WA,[],2) - 1);
% max(WA,[],2)
ISI = ISI + sum(sum(WA,1)./max(WA,[],1) - 1);
% max(WA,[],1)

ISI = ISI/(2*N*(N-1));

end