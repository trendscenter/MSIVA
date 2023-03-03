function combinatorial_optim(O, varargin)
% Solve combinatorial optimization

if ~isempty(varargin)
%     w0 = O.greedysearch_iva(varargin{1});
%     w0 = O.greedysearch_flex(varargin{1});
    w0 = O.greedysearch(varargin{1});
else
%     w0 = O.greedysearch_iva();
%     w0 = O.greedysearch_flex();
    w0 = O.greedysearch();
end

% w0 = O.sub_perm_analysis(w0);

ut = utils;
stdY1 = std(O.Y{1},[],2);
stdY2 = std(O.Y{2},[],2);
w0_us = ut.unstackW(w0,O.M,O.C,O.V);
w0_ = ut.stackW({diag(pi/sqrt(3)./stdY1)*w0_us{1},diag(pi/sqrt(3)./stdY2)*w0_us{2}});

O.objective(w0_);

end