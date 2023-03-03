close all; clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCV info:
tag = 'ISA';
cas = num2str(14);
ISAcaselist = {'5','6','7','9','10','11','13'};
% Features (mixing matrix) info:
casA = 'fake';
r = 2;
%r = ${suf};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load :
f = filesep;
load(fullfile(['..' f 'data' f 'Sgt' f 'case' cas], ...
    ['SCV_' tag '_case' cas '_r' num2str(r,'%03d') '.mat']))

% Generate Features (mixing matrix):
rng(1981+r)
my_cond = 3;
tagA = ['sqcond' num2str(my_cond)];
A = cell(1,max(M));
for mm = M
    [u,s,v] = svd(2*rand(size(Sgt{1},1)+0, size(Sgt{1},1))-1, 'econ');
    A{mm} = u*diag(diag(s)+((max(diag(s))-my_cond*min(diag(s)))/(my_cond-1)))*v';
%     A{ff} = u*v';
end
% for ff = M
%     cond(A{ff})
% end
% A{1} = A{2}(:,1:4);
% A = A(1);

% Generate Mixture
for mm = M
    d_ = sum([S{M}],2);
    dff = diag(d_)*S{mm};
    D = diag(1./sqrt(sum(dff(d_ > 0,:))+1));
    X{mm} = A{mm}*(D*Sgt{mm});
end

% Ground-truth W (cell format)
Wgt = cell(size(A));
for mm = M
    [u,s,v] = svd(A{mm},'econ');
    Wgt{mm} = v*diag(1./diag(s))*u';
end

% Ground-truth W (vector format)
wgt = stackW(Wgt);

total = 10; % Numbere of random starts

%% Run with multiple restarts
delta_f = 1e-5;
max_try = 1;
tp = 'con_bound_';

w0 = [];
W0 = cell(size(A));
for mm = M
    W0{mm} = zeros(size(A{mm}'));
end
w0 = stackW(W0);

myMISA = MISA(w0, M, S, X);

rng(1980 + r)
scurr = rng;
for rr = 1:total
    
    rng(scurr)
    % Random starting pointw
%     l = sqrt(3./V);
    for mm = M
        [u, s, v] = svd(randn(size(A{mm}')),'econ');
        W0{mm} = u*v';
%         W0{ff} = l(ff)*2*(rand(size(A{ff}')) - .5);
    end
    scurr = rng;
    w0 = stackW(W0);
    
    wout = w0;
    Wout = W0;
    
    ll = 1;
    mISIlist = cell(ll);
    fvallist = cell(ll);
    ef = cell(ll);
    outlist = cell(ll);
    
    myMISAobjective = @(x) myMISA.objective_(x);
    
    for ct = 2
        
        mISIlist{ll} = [mISIlist{ll} MISI(Wout,A,myMISA.S)];
        fvallist{ll} = [fvallist{ll} myMISAobjective(wout)];
        ef{ll} = [ef{ll} NaN];
        outlist{ll} = [outlist{ll} struct([])];
        
        optprob = getoptprob4(wout, myMISAobjective);
        [wout,fval,exitflag,output] = fmincon(optprob);
        Wout = myMISA.unstackW(wout,myMISA.M,myMISA.C,myMISA.V);
        
        mISIlist{ll} = [mISIlist{ll} MISI(Wout,A,myMISA.S)];
        fvallist{ll} = [fvallist{ll} myMISAobjective(wout)];
        ef{ll} = [ef{ll} exitflag];
        outlist{ll} = [outlist{ll} output];
        
        if length(myMISA.M) ~= 1 || ... % More than 1 dataset
           myMISA.K ~= myMISA.C(myMISA.M(1)) || ... %
           prod(diag(myMISA.S{myMISA.M(1)}(myMISA.nes,myMISA.C(myMISA.M(1)))) == ones(myMISA.C(myMISA.M(1)),1)) ~= 1
            
            woutT = wout;
            WoutT = Wout;
            exitflagT = exitflag;
            outputT = output;
            
            S_ = myMISA.S;
            for mm = myMISA.M
                S_{mm} = eye(myMISA.C(mm));
                woutT = myMISA.stackW(WoutT(mm));
                myMISA.update(S_,mm);
                
                optprob = getoptprob4(woutT, myMISAobjective);
                % optprob.options.MaxIter = 2;
                [woutT,fvalT,exitflagT,outputT] = fmincon(optprob);
                W = myMISA.unstackW(woutT,myMISA.M,myMISA.C,myMISA.V);
                WoutT{mm} = W{mm};
            end
            
            woutT = myMISA.stackW(WoutT);
            myMISA.update(S,M);
            
            mISIlist{ll} = [mISIlist{ll} MISI(WoutT,A,myMISA.S)];
            fvallist{ll} = [fvallist{ll} myMISAobjective(woutT)];
            ef{ll} = [ef{ll} exitflag];
            outlist{ll} = [outlist{ll} output];
            
            S_ = myMISA.S;
            for mm = myMISA.M
                S_{mm} = eye(myMISA.C(mm));
                woutT = myMISA.stackW(WoutT(mm));
                myMISA.update(S_,mm);
                
                for cc = 1:size(S_{mm}, 2)%randperm(size(S{ff}, 2))
                    current = find(S_{mm}(:,cc));
                    S_{mm}(:,cc) = zeros(size(S_{mm},1),1);
                    misa_values = [];
                    for ss = 1:size(S_{mm},1)
                        if ss ~= 1
                            S_{mm}(ss-1,cc) = 0;
                        end
                        S_{mm}(ss,cc) = 1;
                        myMISA.update(S_,mm);
                        misa_values(ss) = myMISAobjective(woutT);
                    end
                    S_{mm}(ss,cc) = 0;
                    S_{mm} = [S_{mm}; zeros(1,size(S_{mm},2))];
                    S_{mm}(end,cc) = 1;
                    myMISA.update(S_,mm);
                    misa_values(ss+1) = myMISAobjective(woutT);
                    [~,ix] = min(misa_values);
                    if ix ~= current && abs(diff(misa_values([ix current]))) < 1e-8
                        ix = current;
                    end
                    if ix < (ss+1)
                        S_{mm} = S_{mm}(1:(end-1),:);
                        S_{mm}(ix,cc) = 1;
                    end
                    S_{mm} = S_{mm}(sum(S_{mm},2) ~= 0,:);
                end
                shuff = match_subspaces(S{mm},S_{mm});
                S_{mm} = S_{mm}(:,shuff);
                myMISA.update(S_,mm);
                
                WoutT{mm} = WoutT{mm}(shuff,:);
                Wout{mm} = Wout{mm}(shuff,:);
                W0{mm} = W0{mm}(shuff,:);
            end
            
            myMISA.update(S,M);
            woutT = myMISA.stackW(WoutT);
            wout = myMISA.stackW(Wout);
            w0 = myMISA.stackW(W0);
            
            mISIlist{ll} = [mISIlist{ll} MISI(Wout,A,myMISA.S)];
            fvallist{ll} = [fvallist{ll} myMISAobjective(wout)];
            ef{ll} = [ef{ll} NaN];
            outlist{ll} = [outlist{ll} struct([])];
            
            subspace_perms = gen_list_Gperms(S,myMISA.M);
            spval = [];
            S_ = myMISA.S;
            for kk = 1:size(subspace_perms{myMISA.M(1)},1)
                for mm = myMISA.M
                    S_{mm} = S{mm}(subspace_perms{mm}(kk,:),:);
                end
                myMISA.update(S_,M);
                spval(kk) = myMISAobjective(wout);
            end
            [~, the_min] = min(spval);
            for mm = myMISA.M
                S_{mm} = S{mm}(subspace_perms{mm}(the_min,:),:);
            end
            myMISA.update(S_,M);
            
            shuff = cell(1,max(myMISA.M));
            for mm = myMISA.M
                shuff{mm} = [];
                for ss = 1:size(myMISA.S{myMISA.M(1)},1)
                    shuff{mm} = [shuff{mm} find(myMISA.S{mm}(ss,:))];
                end
            end
            myMISA.update(S,M);
            
            for mm = M
                Wout{mm} = Wout{mm}(shuff{mm},:);
                W0{mm} = W0{mm}(shuff{mm},:);
            end
            wout = stackW(Wout);
            w0 = stackW(W0);
            
        end
        
    end
    
    mISIlist{ll} = [mISIlist{ll} MISI(Wout,A,myMISA.S)];
    fvallist{ll} = [fvallist{ll} myMISAobjective(wout)];
    ef{ll} = [ef{ll} NaN];
    outlist{ll} = [outlist{ll} struct([])];
    
    optprob = getoptprob4(wout, myMISAobjective);
    [wout,fval,exitflag,output] = fmincon(optprob);
    Wout = myMISA.unstackW(wout,myMISA.M,myMISA.C,myMISA.V);
    
    mISIlist{ll} = [mISIlist{ll} MISI(Wout,A,myMISA.S)];
    fvallist{ll} = [fvallist{ll} myMISAobjective(wout)];
    ef{ll} = [ef{ll} exitflag];
    outlist{ll} = [outlist{ll} output];
    
%     MISA_ICA_perm;
%     MISA_ICA_perm;
    
    mISI = MISI(Wout,A,myMISA.S);
    save(['..' f 'output' f 'S_case' cas '_A_' casA f 'r' num2str(r,'%03d') ...
            f tp f tag '_' tagA '_' num2str(rr,'%03d') '.mat'], ...
            'scurr', 'w0', 'wout', 'mISI', 'fval', 'exitflag', 'output', ...
            'mISIlist', 'fvallist', 'ef', 'outlist')
end

%% ISA (DLahat)
tp = 'dlahat_';
rng(1980 + r)
scurr = rng;
rp = 100;
B = 5;
spb = floor(N/B)-1;
bs = 1:spb:(N-1);
lbs = length(bs)-1;

for rr = 1 % NO POINT RUNNING MULTIPLE STARTS: ALWAYS STARTS FROM IDENTITY MATRIX
    
    Rxx = cell(1,rp*lbs);
    
    for pp = 1:rp
        seq = randperm(N);
        for bb = 1:lbs
            if bb < lbs
                Rxx{(pp-1)*lbs+bb} = cov(Xdat{M}(:,seq(bs(bb):(bs(bb+1)-1)))');
            else
                Rxx{(pp-1)*lbs+bb} = cov(Xdat{M}(:,seq(bs(bb):end))');
            end
        end
    end
    
    W0 = cell(size(A));
    [W0{M}, ~, rc, phi, swp] = jbd([Rxx{:}], ones(1,rp*lbs), 1e-9, 512, sum([S{:}],2));
    w0 = stackW(W0);
    wout = w0;
    Wout = W0;
    
    fval = MISA(wout);
    mISI = MISI(Wout,A,S);
    save(['..' f 'output' f 'S_case' cas '_A_' casA f 'r' num2str(r,'%03d') ...
            f tp f tag '_' tagA '_' num2str(rr,'%03d') '.mat'], ...
            'scurr', 'w0', 'wout', 'mISI', 'fval')
%     figure
% %     for bb = 1:lbs
% %         imr = abs(Wout{M}*Rxx{bb}*Wout{M}');
%         imr = abs(cov((Wout{M}*Xdat{M})'));
%         imr = (sqrt(diag(diag(imr)))\imr)/sqrt(diag(diag(imr)));
% %         imr = corr(Sgt{M}',(Wout{M}*Xdat{M})');
%         imagesc(imr,max(imr(:))*[0 1])
%         colormap hot
%         axis equal tight
%         drawnow()
%         pause(.5)
% %     end
    
end

%% ISA (QLe)
tp = 'qle_';
rng(1980 + r)
scurr = rng;

for rr = 1:total
    rng(scurr)
    
    % Random starting pointw
    w0 = [];
    W0 = cell(size(A));
    for mm = M
        [u, s, v] = svd((1/sqrt(N-1))*Xdat{mm}, 'econ');
        wM{mm} = s\u';
        Xdat{mm} = wM{mm}*Xdat{mm};
        [u, s, v] = svd(randn(size(A{mm}')), 'econ');
        W0{mm} = u*v';
    end
    scurr = rng;
    w0 = stackW(W0);
    
    Wout = cell(size(A));
    Wout{M} = isa_est(W0{:}, Xdat{:}, wM{:}, C(M), C(M)/K);
    Wout{M} = Wout{M}*wM{:};
    wout = stackW(Wout);
    fval = MISA(wout);
    mISI = MISI(Wout,A,S);
    
    save(['..' f 'output' f 'S_case' cas '_A_' casA f 'r' num2str(r,'%03d') ...
        f tp f tag '_' tagA '_' num2str(rr,'%03d') '.mat'], ...
        'scurr', 'w0', 'wout', 'mISI', 'fval')
end

%%
mISIlist = [];
fvallist = [];
exitflaglist = [];
group = [];
tp = 'con_bound_';
% tp = 'qle_';
% tp = 'dlahat_';
R = 2;
switch lower(tp)
    case 'qle_'
        Rv = R;
    otherwise
        Rv = 1:R;
end
for r = 1:length(Rv)
    group = [group Rv(r)*ones(1,total)];
    for rr = 1:total
        try
            load(['..' f 'output' f 'S_case' cas '_A_' casA f 'r' num2str(Rv(r),'%03d') ...
                f tp f tag '_' tagA '_' num2str(rr,'%03d') '.mat'], ...
                'mISI', 'fval')
            %             'exitflag')
            mISIlist = [mISIlist mISI];
            fvallist = [fvallist fval];
            %         exitflaglist = [exitflaglist exitflag];
        catch
            mISIlist = [mISIlist 10];
            fvallist = [fvallist 10];
        end
    end
end
mISIlist = 20*log10(mISIlist);

switch lower(tp)
    case 'qle_'
        figure('position',[50   50   300   600])
    otherwise
        figure('position',[50   50   450   600])
end
hold on
B = 15;
owd = -50;
imhist_ = zeros(B,length(Rv));
for r = 1:length(Rv)
%     [c, e] = histme_1D(mISIlist(group == r),[min(mISIlist)-eps max(mISIlist)+eps],B);
    ylim_ = min(owd, floor(min(mISIlist(:))/10)*10);
    [c, e] = histme_1D(mISIlist(group == Rv(r)),[ylim_ 0],B);
    imhist_(:,r) = c;
end
es = e(1:B) + (e(2)-e(1))/2;
imagesc([min(Rv) max(Rv)], [es(1) es(B)],imhist_,[0,total])
axis tight
axis xy
ht = hot;
colormap(ht(64:-1:17,:))
set(gca,'xtick',Rv)
set(gca,'ylim',[e(1) e(end)])
yticks = 0:-5:ylim_;
set(gca,'ytick',yticks(end:-1:1))
% set(gca,'ytick',e)
% set(gca,'yticklabel',cellstr(num2str(e,'%.1f')))
switch lower(tp)
    case 'con_bound_'
        title('MISA (with SOS) solving an ISA problem')
    case 'dlahat_'
        title('2^{nd}-Order MICA solving an ISA problem')
    case 'qle_'
        title('isa\_est() solving an ISA problem')
end
ch = colorbar;
ylabel(ch, 'Counts')
switch lower(cas)
    case '13'
        xlabel('Datasets without correlation')
    case '14'
        xlabel('Datasets with correlation')
end
ylabel('joint ISI [dB]')

hold on
for r = 1:length(Rv)
    plot(group(group == Rv(r)) + .6*rand(1,sum(group==Rv(r)))-.3,mISIlist(group==Rv(r)),'.k')
end

