function T = simISA(Run, part_a, part_b, inn, iaa, ir0, er0)
%close all; clear all; clc
format compact

ut = utils;

fname = ['ISAsim_w0_results_rr' num2str(Run,'%03d') '_ir0' num2str(ir0,'%03d') '_er0' num2str(er0,'%03d') '.csv'];
if exist(fullfile(pwd, 'results', 'sim', 'ISA', fname),'file')
    T = readtable(fullfile(pwd, 'results', 'sim', 'ISA', fname));
    ir0 = T.R0(end);
else
    T = cell2table(cell(0,12));
    T.Properties.VariableNames = {'SNRdB','Acond','Run','R0','Initial',...
        'Alg','Iter','WallTime','Cost','MISI','MMD','MMSE'};
end

% C = 75; V = 8000;
C = [1:5, 5:-1:1, 1:5, 2, 2, 2]; V = 8000;
% S_      = [{[1 2 3];   [4]; [5 6 7];       [8]; [9]}, ...
%            {  [1 2]; [3 4];      []; [5 6 7 8]; [9]}, ...
%            {  [1 2];   [3]; [4 5 6];   [7 8 9];  []}];
% V = [2000 1000 1500];
N = 5250;
Acond = [1 3 7 15];
SNR = [999   9   1   1   1; % source power
         1   1   1   9 999];% noise power
SNR = sum(SNR)./SNR(2,:); %(1+99)/99;

if part_a & part_b
    inna = inn; innb = 1;
    iaaa = iaa; iaab = 1;
    ir0a = ir0; ir0b = ir0;
    er0a = er0; er0b = er0;
else
    inna = inn; innb = inn;
    iaaa = iaa; iaab = iaa;
    ir0a = ir0; ir0b = ir0;
    er0a = er0; er0b = er0;
end

rng(100)
r0seedsa = floor(10000*rand(1,1000));
r0seedsb = floor(10000*rand(1,1000));

if part_a
    nn = 3;
    for aa = iaaa:size(Acond,2)
        rr = Run;
        seed = 1980+rr;
        sim1 = sim_basic_ISA(seed,C,V,N,Acond(aa),SNR(nn));
        seed = rng;
        
        beta = .5; gradtype = 'relative'; sc = true; RElambda = 0;
        data1 = setup_basic_MISAKRE(sim1, beta, gradtype, sc, RElambda);
        
        for r0 = ir0a:er0a%1%10
            fprintf('SNR: %.3f, Acond: %d, Run: %d, R0: %d(%d-%d)\n', ...
                SNR(nn), Acond(aa), rr, r0, ir0b, er0b)
            
            W0 = cell(size(sim1.A));
            for mm = sim1.M
                %                 W0{mm} = randn(size(sim1.A{mm}'));
                rng(r0seedsa(r0))
                [u, s, v] = svd(randn(size(sim1.A{mm}')),'econ');
                W0{mm} = u*v';
            end
            
            %             Mold = data1.M;
            %             Sold = data1.S;
            %             bold = data1.beta;
            %             lold = data1.lambda;
            %             eold = data1.eta;
            %             old_d_k = data1.d_k;
            %             oldsc = data1.sc; % Save scale-control parameter
            %             S_ = Sold;
            %
            %         for mm = Mold
            %             % Set S to sparse identity matrix
            %             S_{mm} = sparse(1:data1.C(mm),1:data1.C(mm),ones(data1.C(mm),1), ...
            %                 data1.C(mm),data1.C(mm),data1.C(mm)); %eye(O.C(mm));
            %             % redefine parameters: accounts for changing subspace structure/size
            %             ixe = cumsum(cell2mat(old_d_k(mm)));
            %             ixb = [1; (ixe(1:end-1)+1)];
            %             b = []; l = []; e = [];
            %             for kk = 1:length(ixe)
            %                 b(ixb(kk):ixe(kk)) = bold(1);
            %                 l(ixb(kk):ixe(kk)) = lold(1);
            %                 e(ixb(kk):ixe(kk)) = eold(1);
            %             end
            %             data1.update(S_,mm,b',l',e'); % Update S and set M = mm
            w0 = ut.stackW(W0(mm));
            optprob = ut.getop(w0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 1e-7);
            optprob.options.MaxIter = 1;
%             data1.setREapproach('WT')
            [woutW0,fval,exitflag,output] = fmincon(optprob);
            optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 10^floor(log10(output.firstorderopt))/80);
            %         optprob.options.MaxFunEvals = 20;
            %         [woutW0,fval,exitflag,output] = fmincon(optprob);
            %         optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, optprob.options.TolX);
            % %         optprob.options.MaxFunEvals = 15;
%             data1.setREapproach('PINV')
            [woutW0,fval,exitflag,output] = fmincon(optprob);
            %         end
            %         data1.update(Sold,Mold,bold,lold,eold);
            
            tic
            data2 = MISAKRE(ut.stackW(cellfun(@(c) eye(c), num2cell(data1.C), 'Un', 0)), ...
                data1.M, data1.S, data1.Y, data1.beta, data1.eta, data1.lambda, data1.gradtype, data1.sc, data1.preX, ...
                data1.REtype, data1.REapproach, data1.RElambda, ...
                data1.REref, data1.REreftype, data1.REreflambda, data1.rC);
            
            barrier = 1; m = 15;
            optprob = ut.getop(ut.stackW(data2.W), @(x) data2.objective(x), [], barrier, {'lbfgs' m});
            [wout,fval,exitflag,output] = fmincon(optprob);
            % [w, fval, iters] = myline_search_batch(ut.stackW(data2.W),data2);
            t = toc;
            
            %         figure, imagesc(corr(data2.Y{3}',sim1.Y{3}'), [-1 1])
            %         figure, imagesc(corr(data2.Y{2}',sim1.Y{2}'), [-1 1])
            %         figure, imagesc(corr(data2.Y{1}',sim1.Y{1}'), [-1 1])
            %         ut.MISI(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0), sim1.A, data1.S)
            
            data2.updatesc(false)
            data2.objective(ut.stackW(data2.W));
            data2.updatesc(true)
            data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
            
            T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
                {'MISA w/o GP'}, output.iterations, t, data2.objective_(), ...
                data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
            T2.Properties.VariableNames = T.Properties.VariableNames;
            T = [T;T2];
            
            fprintf('MISA w/o GP: %.5f\n', data1.MISI(sim1.A))
            
            data1.objective(woutW0); %reset
            aux = {data2.W; data2.objective(ut.stackW(data2.W))};
            for ct = 2:3
                data2.combinatorial_optim()
                optprob = ut.getop(ut.stackW(data2.W), @(x) data2.objective(x), [], barrier, {'lbfgs' m});
                [wout,fval,exitflag,output] = fmincon(optprob);
                aux(:,ct) = {data2.W; data2.objective_()};
            end
            [~, ix] = min([aux{2,:}]);
            t = toc;
            
            data2.updatesc(false)
            data2.objective(ut.stackW(aux{1,ix}));
            data2.updatesc(true)
            data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
            
            T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
                {'MISA w/ GP'}, output.iterations, t, data2.objective_(), ...
                data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
            T2.Properties.VariableNames = T.Properties.VariableNames;
            T = [T;T2];
            
            fprintf('MISA w/ GP: %.5f\n', data1.MISI(sim1.A))
            
            % JBD - DLahat
            data1.objective(woutW0); %reset
            tic
            data2 = MISAKRE(ut.stackW(cellfun(@(c) eye(c), num2cell(data1.C), 'Un', 0)), ...
                data1.M, data1.S, data1.Y, data1.beta, data1.eta, data1.lambda, data1.gradtype, data1.sc, data1.preX, ...
                data1.REtype, data1.REapproach, data1.RElambda, ...
                data1.REref, data1.REreftype, data1.REreflambda, data1.rC);
            
%             Rxx = cell(1,N-999);
%             seq = randperm(N);
%             for pp = 0:(N-1000)
%                 Rxx{pp+1} = (1/(N-1))*(data2.X{1}(:,seq(1:(end-pp))) * data2.X{1}(:,seq((1+pp):end))');
%             end
            rp = 100;
            B = 5;
            spb = floor(N/B)-1;
            bs = 1:spb:(N-1);
            lbs = length(bs)-1;
            Rxx = cell(1,rp*lbs);
            for pp = 1:rp
                seq = randperm(N);
                for bb = 1:lbs
                    if bb < lbs
                        Rxx{(pp-1)*lbs+bb} = cov(data2.X{1}(:,seq(bs(bb):(bs(bb+1)-1)))');
                    else
                        Rxx{(pp-1)*lbs+bb} = cov(data2.X{1}(:,seq(bs(bb):end))');
                    end
                end
            end
            %[W, ~, rc, phi, swp] = jbd([Rxx{:}], ones(1,N-999), 1e-9, 512, full(sum([data2.S{:}],2)));
            [W, ~, rc, phi, swp] = jbd([Rxx{:}], ones(1,length(Rxx)), 1e-6, 512, full(sum([data2.S{:}],2)));
            t = toc;
            
            data2.updatesc(false)
            data2.objective(ut.stackW({W}));
            data2.updatesc(true)
            data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
            
            T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
                {'JBD SOS'}, output.iterations, t, data2.objective_(), ...
                data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
            T2.Properties.VariableNames = T.Properties.VariableNames;
            T = [T;T2];
            
            fprintf('JBD SOS: %.5f\n', data1.MISI(sim1.A))
            
            writetable(T, fullfile(pwd, 'results', 'sim', 'ISA', fname))
            
        end
%         ir0a = 1;
    end
    iaaa = 1;
end

if part_b
    aa = 3;
    for nn = innb:size(SNR,2)
        rr = Run;
        seed = 1990+rr;
        sim1 = sim_basic_ISA(seed,C,V,N,Acond(aa),SNR(nn));
        seed = rng;
        
        beta = .5; gradtype = 'relative'; sc = true; RElambda = 0;
        data1 = setup_basic_MISAKRE(sim1, beta, gradtype, sc, RElambda);
        
        for r0 = ir0b:er0b%1%10
            fprintf('SNR: %.3f, Acond: %d, Run: %d, R0: %d(%d-%d)\n', ...
                SNR(nn), Acond(aa), rr, r0, ir0b, er0b)
            
            W0 = cell(size(sim1.A));
            for mm = sim1.M
                %                 W0{mm} = randn(size(sim1.A{mm}'));
                rng(r0seedsb(r0))
                [u, s, v] = svd(randn(size(sim1.A{mm}')),'econ');
                W0{mm} = u*v';
            end
            
            %             Mold = data1.M;
            %             Sold = data1.S;
            %             bold = data1.beta;
            %             lold = data1.lambda;
            %             eold = data1.eta;
            %             old_d_k = data1.d_k;
            %             oldsc = data1.sc; % Save scale-control parameter
            %             S_ = Sold;
            %
            %         for mm = Mold
            %             % Set S to sparse identity matrix
            %             S_{mm} = sparse(1:data1.C(mm),1:data1.C(mm),ones(data1.C(mm),1), ...
            %                 data1.C(mm),data1.C(mm),data1.C(mm)); %eye(O.C(mm));
            %             % redefine parameters: accounts for changing subspace structure/size
            %             ixe = cumsum(cell2mat(old_d_k(mm)));
            %             ixb = [1; (ixe(1:end-1)+1)];
            %             b = []; l = []; e = [];
            %             for kk = 1:length(ixe)
            %                 b(ixb(kk):ixe(kk)) = bold(1);
            %                 l(ixb(kk):ixe(kk)) = lold(1);
            %                 e(ixb(kk):ixe(kk)) = eold(1);
            %             end
            %             data1.update(S_,mm,b',l',e'); % Update S and set M = mm
            w0 = ut.stackW(W0(mm));
            optprob = ut.getop(w0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 1e-7);
            optprob.options.MaxIter = 1;
%             data1.setREapproach('WT')
            [woutW0,fval,exitflag,output] = fmincon(optprob);
            optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 10^floor(log10(output.firstorderopt))/4);
            %         optprob.options.MaxFunEvals = 20;
            %         [woutW0,fval,exitflag,output] = fmincon(optprob);
            %         optprob = ut.getop(woutW0, @(x) data1.opt_RE(x), [], 0.1, {'lbfgs' 5}, 10^floor(log10(output.firstorderopt))/40);
            % %         optprob.options.MaxFunEvals = 15;
%             data1.setREapproach('PINV')
            [woutW0,fval,exitflag,output] = fmincon(optprob);
            %         end
            %         data1.update(Sold,Mold,bold,lold,eold);
            
            tic
            data2 = MISAKRE(ut.stackW(cellfun(@(c) eye(c), num2cell(data1.C), 'Un', 0)), ...
                data1.M, data1.S, data1.Y, data1.beta, data1.eta, data1.lambda, data1.gradtype, data1.sc, data1.preX, ...
                data1.REtype, data1.REapproach, data1.RElambda, ...
                data1.REref, data1.REreftype, data1.REreflambda, data1.rC);
            
            barrier = 1; m = 15;
            optprob = ut.getop(ut.stackW(data2.W), @(x) data2.objective(x), [], barrier, {'lbfgs' m});
            [wout,fval,exitflag,output] = fmincon(optprob);
            % [w, fval, iters] = myline_search_batch(ut.stackW(data2.W),data2);
            t = toc;
            
            %         figure, imagesc(corr(data2.Y{3}',sim1.Y{3}'), [-1 1])
            %         figure, imagesc(corr(data2.Y{2}',sim1.Y{2}'), [-1 1])
            %         figure, imagesc(corr(data2.Y{1}',sim1.Y{1}'), [-1 1])
            %         ut.MISI(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0), sim1.A, data1.S)
            
            data2.updatesc(false)
            data2.objective(ut.stackW(data2.W));
            data2.updatesc(true)
            data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
            
            T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
                {'MISA w/o GP'}, output.iterations, t, data2.objective_(), ...
                data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
            T2.Properties.VariableNames = T.Properties.VariableNames;
            T = [T;T2];
            
            fprintf('MISA w/o GP: %.5f\n', data1.MISI(sim1.A))
            
            data1.objective(woutW0); % reset
            aux = {data2.W; data2.objective(ut.stackW(data2.W))};
            for ct = 2:3
                data2.combinatorial_optim()
                optprob = ut.getop(ut.stackW(data2.W), @(x) data2.objective(x), [], barrier, {'lbfgs' m});
                [wout,fval,exitflag,output] = fmincon(optprob);
                aux(:,ct) = {data2.W; data2.objective_()};
            end
            [~, ix] = min([aux{2,:}]);
            t = toc;
            
            data2.updatesc(false)
            data2.objective(ut.stackW(aux{1,ix}));
            data2.updatesc(true)
            data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
            
            T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
                {'MISA w/ GP'}, 0, t, data2.objective_(), ...
                data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
            T2.Properties.VariableNames = T.Properties.VariableNames;
            T = [T;T2];
            
            fprintf('MISA w/ GP: %.5f\n', data1.MISI(sim1.A))
            
            % JBD - DLahat
            data1.objective(woutW0); %reset
            tic
            data2 = MISAKRE(ut.stackW(cellfun(@(c) eye(c), num2cell(data1.C), 'Un', 0)), ...
                data1.M, data1.S, data1.Y, data1.beta, data1.eta, data1.lambda, data1.gradtype, data1.sc, data1.preX, ...
                data1.REtype, data1.REapproach, data1.RElambda, ...
                data1.REref, data1.REreftype, data1.REreflambda, data1.rC);
            
%             Rxx = cell(1,N-999);
%             seq = randperm(N);
%             for pp = 0:(N-1000)
%                 Rxx{pp+1} = (1/(N-1))*(data2.X{1}(:,seq(1:(end-pp))) * data2.X{1}(:,seq((1+pp):end))');
%             end
            rp = 100;
            B = 5;
            spb = floor(N/B)-1;
            bs = 1:spb:(N-1);
            lbs = length(bs)-1;
            Rxx = cell(1,rp*lbs);
            for pp = 1:rp
                seq = randperm(N);
                for bb = 1:lbs
                    if bb < lbs
                        Rxx{(pp-1)*lbs+bb} = cov(data2.X{1}(:,seq(bs(bb):(bs(bb+1)-1)))');
                    else
                        Rxx{(pp-1)*lbs+bb} = cov(data2.X{1}(:,seq(bs(bb):end))');
                    end
                end
            end
            %[W, ~, rc, phi, swp] = jbd([Rxx{:}], ones(1,N-999), 1e-9, 512, full(sum([data2.S{:}],2)));
            [W, ~, rc, phi, swp] = jbd([Rxx{:}], ones(1,length(Rxx)), 1e-6, 512, full(sum([data2.S{:}],2)));
            t = toc;
            
            data2.updatesc(false)
            data2.objective(ut.stackW({W}));
            data2.updatesc(true)
            data1.objective(ut.stackW(cellfun(@(d2w,d1w) d2w*d1w,data2.W,data1.W,'Un',0)));
            
            T2 = cell2table({10*log10(SNR(nn)), Acond(aa), rr, r0, {'W0'}, ...
                {'JBD SOS'}, 0, t, data2.objective_(), ...
                data1.MISI(sim1.A), data1.MMD(sim1.A), data1.MMSE(sim1.Y) });
            T2.Properties.VariableNames = T.Properties.VariableNames;
            T = [T;T2];
            
            fprintf('JBD SOS: %.5f\n', data1.MISI(sim1.A))
            
            writetable(T, fullfile(pwd, 'results', 'sim', 'ISA', fname))
            
        end
%         ir0b = 1;
    end
    innb = 1;
end

