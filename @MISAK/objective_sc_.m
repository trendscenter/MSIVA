function [J, gJ] = objective_sc_(O)

% Compute cost function: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MISA cost/objective function:
JE = 0; % 1st TERM, a
JF = 0; % 1st TERM, b
JC = 0; % 2nd TERM
JD = 0; % 3rd TERM
fc = 0; % Constants
if nargout > 1
    % Initialize MISA cost/objective gradient:
    % 1st TERM:
    gJ = cell(1,max(O.M));
    switch O.gradtype
        case 'regular'
            gJ(O.M) = cellfun(@(w) zeros(size(w)), O.W(O.M),'Un',0);
        case 'relative'
            gJ(O.M) = cellfun(@(w) zeros(size(w,1)), O.W(O.M),'Un',0);
        otherwise
            gJ(O.M) = cellfun(@(w) zeros(size(w)), O.W(O.M),'Un',0);
    end
end

% 1st and 2nd TERMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sum of subspace entropies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t1 = cell(find(O.nes', 1, 'last' ),1);
for kk = find(O.nes')
    y_sub = zeros(O.d(kk),O.N);
    tot = 0;
    for mm = O.M
        ix = tot + (1:O.d_k{mm}(kk));
        if ~isempty(ix)
            y_sub(ix,:) = O.Y{mm}(logical(O.S{mm}(kk,:)),:);
        end
        tot = tot + O.d_k{mm}(kk);
    end
%     y_sub = cell2mat(cellfun(@(y,s) y(logical(s(kk,:)'),:), O.Y(O.M), O.S(O.M),'Un',0)');
    yyT = y_sub*y_sub';
    g_k = diag(yyT).^(-.5);
    g2_k = diag(yyT).^(-1);
    g_kInv = diag(yyT).^(.5);
    ybar_sub = bsxfun(@times, g_kInv, y_sub);
%     A = pinv(y_sub)'; This is insanely slow...
    % The following 3 lines are a more efficient way to compute pinv():
%     [t1, D] = eig(y_sub*y_sub');
%     D = diag(D);
%     A = (bsxfun(@times,t1,1./D')*t1')*y_sub;
    % Unfortunatelly, scale control requires to save the inverse of y_sub*y_sub':
%     yyTInv = bsxfun(@times,t1,1./D')*t1';
    yyTInv = yyT\eye(O.d(kk));
    %clear yyT
    A = yyTInv * ybar_sub;
    z_k = sum(ybar_sub .* A, 1);
    z_k_beta = z_k.^O.beta(kk);
    
    JE = JE + O.lambda(kk) * mean(z_k_beta);
    if O.eta(kk) ~= 1
        JF = JF + (1 - O.eta(kk)) * mean(log(z_k));
    end
    JC = JC + sum(log(eig(bsxfun(@times,g_k,bsxfun(@times,yyT,g_k)))));
    % Gradient:
    if nargout > 1
        % Number of dimensions in each modality (for subspace kk):
        z_k = ((2*O.beta(kk)*O.lambda(kk)/O.N)*z_k_beta + (2*(1-O.eta(kk))/O.N))./z_k;
        B = bsxfun(@times, A, z_k);
        C = -(B * A');
        cc = sub2ind(size(C),1:length(C),1:length(C))';
        C(cc) = C(cc) + g2_k.*sum(B.*ybar_sub,2);
        yyTInv(cc) = yyTInv(cc) - g2_k;
        C = C + yyTInv;
        B = bsxfun(@times, g_kInv, B);
        
%         d_k = cell2mat(cellfun(@(s) sum(s(kk,:)), O.S(O.M),'Un',0));
        t1 = cell(size(O.W));
        switch O.gradtype
            case 'regular'
                tot = 0;
                for mm = O.M
                    ix = tot + (1:O.d_k{mm}(kk));
                    if ~isempty(ix)
                        t1{mm} = (B(ix,:) + C(ix,:)*y_sub) * O.X{mm}';
%                         t1{kk}{mm} = (B(ix,:) + C(ix,:)*y_sub) * O.X{mm}';
                    else
                        t1{mm} = zeros(0,O.V(mm));
                    end
                    tot = tot + O.d_k{mm}(kk);
                end
%                 t1{kk} = cellfun(@(r,x) r*x', mat2cell((B + C*y_sub), d_k)', O.X(O.M),'Un',0);
            case 'relative'
                tot = 0;
                for mm = O.M
                    ix = tot + (1:O.d_k{mm}(kk));
                    if ~isempty(ix)
                        t1{mm} = (B(ix,:) + C(ix,:)*y_sub) * O.Y{mm}';
%                         t1{kk}{mm} = (B(ix,:) + C(ix,:)*y_sub) * O.Y{mm}';
                    else
                        t1{mm} = zeros(0,O.C(mm));
                    end
                    tot = tot + O.d_k{mm}(kk);
                end
%                 t1{kk} = cellfun(@(r,y) r*y', mat2cell((B + C*y_sub), d_k)', O.Y(O.M),'Un',0);
            otherwise
                tot = 0;
                for mm = O.M
                    ix = tot + (1:O.d_k{mm}(kk));
                    if ~isempty(ix)
                        t1{mm} = (B(ix,:) + C(ix,:)*y_sub) * O.X{mm}';
%                         t1{kk}{mm} = (B(ix,:) + C(ix,:)*y_sub) * O.X{mm}';
                    else
                        t1{mm} = zeros(0,O.V(mm));
                    end
                    tot = tot + O.d_k{mm}(kk);
                end
%                 t1{kk} = cellfun(@(r,x) r*x', mat2cell((B + C*y_sub), d_k)', O.X(O.M),'Un',0);
        end
        for mm = O.M
            switch O.gradtype
                case 'regular'
                    gJ{mm}(logical(O.S{mm}(kk,:)'),:) = ...
                        gJ{mm}(logical(O.S{mm}(kk,:)'),:) + t1{mm};
                case 'relative'
                    gJ{mm}(logical(O.S{mm}(kk,:)'),:) = ...
                        gJ{mm}(logical(O.S{mm}(kk,:)'),:) + t1{mm};
                otherwise
                    gJ{mm}(logical(O.S{mm}(kk,:)'),:) = ...
                        gJ{mm}(logical(O.S{mm}(kk,:)'),:) + t1{mm};
            end
        end
    end
end
% if nargout > 1
%     for kk = find(O.nes')
%         for mm = O.M
%             switch O.gradtype
%                 case 'regular'
%                     gJ{mm}(logical(O.S{mm}(kk,:)'),:) = ...
%                         gJ{mm}(logical(O.S{mm}(kk,:)'),:) + t1{kk}{mm};
%                 case 'relative'
%                     gJ{mm}(logical(O.S{mm}(kk,:)'),:) = ...
%                         gJ{mm}(logical(O.S{mm}(kk,:)'),:) + t1{kk}{mm};
%                 otherwise
%                     gJ{mm}(logical(O.S{mm}(kk,:)'),:) = ...
%                         gJ{mm}(logical(O.S{mm}(kk,:)'),:) + t1{kk}{mm};
%             end
%         end
%     end
% end

JC = JC/2;

%clear A B C t1 z_k y_sub %cc d_k D

% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fc = .5*log(pi)*sum(O.d) + sum(gammaln(O.nu)) - ...
    sum(gammaln(.5*O.d)) - sum(O.nu.*log(O.lambda)) - sum(log(O.beta));

% 3rd TERM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proxy of log abs det W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    for mm = O.M
        [rr, cc] = size(O.W{mm});
        if rr == cc
            JD = JD - log(abs(det(O.W{mm})));
            switch O.gradtype
                case 'regular'
                    gJ{mm} = gJ{mm} - inv(O.W{mm})';
                case 'relative'
                    gJ{mm} = (gJ{mm} - eye(size(gJ{mm})))*O.W{mm};
                otherwise
                    gJ{mm} = gJ{mm} - inv(O.W{mm})';
            end
        else
            [t1, D] = eig(O.W{mm}*O.W{mm}');
            D = diag(D);
            JD = JD - sum(log(abs(sqrt(D))));
            switch O.gradtype
                case 'regular'
                    gJ{mm} = gJ{mm} - (bsxfun(@times,t1,1./D')*t1')*O.W{mm};
                case 'relative'
                    gJ{mm} = (gJ{mm} - eye(size(gJ{mm})))*O.W{mm};
                otherwise
                    gJ{mm} = gJ{mm} - (bsxfun(@times,t1,1./D')*t1')*O.W{mm};
            end
%             clear t1 D
        end
    end
else
    for mm = O.M
        [rr, cc] = size(O.W{mm});
        if rr == cc
            JD = JD - log(abs(det(O.W{mm})));
        else
            D = eig(O.W{mm}*O.W{mm}');
            JD = JD - sum(log(abs(sqrt(D))));
%             clear D
        end
    end
end

J = JE + JF + JC + JD + fc;
% J=-J;
% for mm = O.M
%     gJ{mm} = -gJ{mm};
% end
