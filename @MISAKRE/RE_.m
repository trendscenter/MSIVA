function [c, DC] = RE_(O)

c = {};%cell(size(O.W));
% 
if nargout > 1
    DC = {};%cellfun(@(w) zeros(size(w)), O.W(O.M),'Un',0);
end

for mm = O.M
    if ~isempty(O.REtype_)
        switch O.REtype_
            case true                       % 'NMSE' (log of)
                switch O.REapproach_
                    case true                   % 'PINV'
                        %                     if O.C(mm) ~= O.V(mm)
                        A = (O.W{mm}*O.W{mm}')\O.W{mm};
                        z2 = A'*O.Y{mm} - O.X{mm};
                        %     x2 = sum(sum(O.X{mm}.^2));
                        %                         [mean(sum(z2.^2)) O.Xnorm(mm) mean(sum(z2.^2))/O.Xnorm(mm)]
                        c{1,mm} = mean(sum(z2.^2))/O.Xnorm(mm); % default
                        %                         c{mm} = sum(z2(:).^2)/O.Xnorm(mm); % for log scale: log(sum_M Rmn)
                        % the following 2 lines are for log scale: sum_M log(Rmn)
                        %                         sRmn = sum(z2(:).^2);
                        %                         c{mm} = (log(sRmn) - log(O.N) - log(O.Xnorm(mm)))/log(10);
                        
                        if nargout > 1
                            AB = (A*O.X{mm})*z2';
                            AB = (2/(O.N*O.Xnorm(mm)))*(AB + (A*z2)*O.X{mm}'); % default
                            %                             AB = (AB + (A*z2)*O.X{mm}'); % for log scale: log(sum_M Rmn)
                            % the following line is for log scale: sum_M log(Rmn)
                            %                             AB = (2/(log(10)*sRmn))*(AB + (A*z2)*O.X{mm}');
                            DC{1,mm} = AB - (AB*A')*O.W{mm};
                        end
                        %                     end
                        
                    case false                  % 'WT'
                        z2 = O.W{mm}'*O.Y{mm} - O.X{mm};
                        c{1,mm} = mean(sum(z2.^2))/O.Xnorm(mm); % default
                        %                     c{mm} = sum(z2(:).^2)/O.Xnorm(mm); % for log scale: log(sum_M Rmn)
                        % the following 2 lines are for log scale: sum_M log(Rmn)
                        %                     sRmn = sum(z2(:).^2);
                        %                     c{mm} = (log(sRmn) - log(O.N) - log(O.Xnorm(mm)))/log(10);
                        
                        if nargout > 1
                            AB = (O.Y{mm})*z2';
                            DC{1,mm} = (2/(O.N*O.Xnorm(mm)))*(AB + (O.W{mm}*z2)*O.X{mm}'); % default
                            %                         DC{mm} = (AB + (O.W{mm}*z2)*O.X{mm}'); % for log scale: log(sum_M Rmn)
                            % the following line is for log scale: sum_M log(Rmn)
                            %                         DC{mm} = (2/(log(10)*sRmn))*(AB + (O.W{mm}*z2)*O.X{mm}');
                        end
                end
            case false                      % 'MSE'
                switch O.REapproach_
                    case true                   % 'PINV'
                        %                     if O.C(mm) ~= O.V(mm)
                        A = (O.W{mm}*O.W{mm}')\O.W{mm};
                        z2 = A'*O.Y{mm} - O.X{mm};
                        %      z2 = (A'*O.Y{mm} - O.X{mm})./O.X{mm};
                        %                         [mean(sum(z2.^2)) O.Xnorm(mm) mean(sum(z2.^2))/O.Xnorm(mm)]
                        c{1,mm} = mean(sum(z2.^2)); % default
                        %                         c{mm} = sum(z2(:).^2);
                        
                        if nargout > 1
                            %      z2 = z2./O.X{mm};
                            AB = (A*O.X{mm})*z2';
                            AB = (2/O.N)*(AB + (A*z2)*O.X{mm}'); % default
                            %                             AB = (AB + (A*z2)*O.X{mm}');
                            DC{1,mm} = AB - (AB*A')*O.W{mm};
                        end
                        %                     end
                        
                    case false                  % 'WT'
                        z2 = O.W{mm}'*O.Y{mm} - O.X{mm};
                        c{1,mm} = mean(sum(z2.^2)); % default
                        %                     c{mm} = sum(z2(:).^2);
                        
                        if nargout > 1
                            AB = (O.Y{mm})*z2';
                            DC{1,mm} = (2/O.N)*(AB + (O.W{mm}*z2)*O.X{mm}'); % default
                            %                         DC{mm} = (AB + (O.W{mm}*z2)*O.X{mm}');
                        end
                end
        end
    end
    if ~isempty(O.REref)
        switch lower(O.REreftype{mm})
            case 'linearreg'
                if ~isempty(O.REref{mm})
                    B = bsxfun(@minus, O.W{mm}(O.rC{mm},:), mean(O.W{mm}(O.rC{mm},:),2));
                    z7 = B*O.W{mm}(O.rC{mm},:)';
                    A = (-z7)\B;
                    e = repmat(mean(O.REref{mm}),O.V(mm),1) - (A'*B)*O.REref{mm} - O.REref{mm};
                    c{2,mm} = mean(sum(e.^2,2))/O.RErefnorm{mm};
                    if nargout > 1
                        AB = (A*O.REref{mm})*e';
                        AB = (2/(O.V(mm)*O.RErefnorm{mm}))*(AB + (A*e)*O.REref{mm}');
                        B = (AB*A')*O.W{mm}(O.rC{mm},:);
                        DC{2,mm} = zeros(O.C(mm),O.V(mm));
                        DC{2,mm}(O.rC{mm},:) = bsxfun(@minus, mean(AB,2), AB) - bsxfun(@plus, B, mean(B));
                    end
%                 else
%                     DC{2,mm} = zeros(O.C(mm),O.V(mm));
                end
            case 'linearreg_nointercept'
                if ~isempty(O.REref{mm})
                    A = (O.W{mm}*O.W{mm}')\O.W{mm};
                    z2 = A'*(O.W{mm}*O.REref{mm}) - O.REref{mm};
                    c{2,mm} = mean(sum(z2.^2,2))/O.RErefnorm{mm};
                    if nargout > 1
                        AB = (A*O.REref{mm})*z2';
                        AB = (2/(O.V(mm)*O.RErefnorm{mm}))*(AB + (A*z2)*O.REref{mm}');
                        DC{2,mm} = AB - (AB*A')*O.W{mm};
                    end
%                 else
%                     DC{2,mm} = zeros(O.C(mm),O.V(mm));
                end
            case 'logisticls'
                if ~isempty(O.REref{mm})
                    B = bsxfun(@minus, O.W{mm}, mean(O.W{mm},2));
                    z7 = B*O.W{mm}';
                    A = (-z7)\B;
                    e = mean(O.REref{mm}) - (A'*B)*O.REref{mm} - O.REref{mm};
                    c{2,mm} = mean(sum(e.^2,2))/O.RErefnorm{mm};
                    if nargout > 1
                        AB = (A*O.REref{mm})*e';
                        AB = (2/(O.V(mm)*O.RErefnorm{mm}))*(AB + (A*e)*O.REref{mm}');
                        B = (AB*A')*O.W{mm};
                        DC{2,mm} = bsxfun(@minus, mean(AB,2), AB) - bsxfun(@plus, B, mean(B));
                    end
%                 else
%                     DC{2,mm} = zeros(O.C(mm),O.V(mm));
                end
            case 'logisticirls'
                if ~isempty(O.REref{mm})
                    
%                 else
%                     DC{2,mm} = zeros(O.C(mm),O.V(mm));
                end
        end
    end
end

end