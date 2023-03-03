%
close all; clear all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCV info:
tag = 'IVA';
cas = num2str(3);
ISAcaselist = {'5','6','7','9','10','11'};
% Features (mixing matrix) info:
casA = 'fake';
r = 2;
%r = ${suf};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load :
f = filesep;
load(fullfile(['..' f 'data' f 'Sgt' f 'case' cas], ...
    ['SCV_' tag '_case' cas '_r' num2str(r,'%03d') '.mat']))

% for mm = M
%     Sgt{mm} = Sgt{mm}(1:5,:);
%     S{mm} = S{mm}(1:5,1:5);
% end

% Generate Features (mixing matrix):
rng(1981+r)
my_cond = 1;
AR = 2;
if AR == 1
    tagA = ['sqcond' num2str(my_cond)];
else
    tagA = ['rectcond' num2str(my_cond)];
end
A = cell(1,max(M));
for ff = M
    cc = size(Sgt{ff},1);
    [u,s,v] = svd(2*rand(cc + [(AR-1)*cc 0])-1, 'econ');
    if my_cond == 1
        A{ff} = u*v';
    else
        A{ff} = u*diag(diag(s)+((max(diag(s))-my_cond*min(diag(s)))/(my_cond-1)))*v';
    end
        %     A{ff} = u*v';
end
% for ff = M
%     cond(A{ff})
% end
% A{1} = A{2}(:,1:4);
% A = A(1);

% Generate Mixture
for ff = M
    d_ = sum([S{M}],2);
    dff = diag(d_)*S{ff};
    D = diag(1./sqrt(sum(dff(d_ > 0,:))+1));
    X{ff} = A{ff}*(D*Sgt{ff});
end

% Ground-truth W (cell format)
Wgt = cell(size(A));
for ff = M
    [u,s,v] = svd(A{ff},'econ');
    Wgt{ff} = v*diag(1./diag(s))*u';
end

% Ground-truth W (vector format)
wgt = stackW(Wgt);
REtype = 'NMSE';
REapproach = 'PINV';
RElambda = 0.5;
myMISARE = MISARE(wgt, M, S, X, REtype, REapproach, RElambda);

%% 
REtype = 'NMSE'; myMISARE.setREtype(REtype);
REapproach = 'PINV'; myMISARE.setREapproach(REapproach);
RElambda = .0001; myMISARE.updateRElambda(RElambda)
optprob = getoptprob4(myMISARE.stackW(Wout), ... %w0, ... %
    @(x)myMISARE.objective(x), @(x) myMISARE.con_RE(x));%@(x) myMISARE.reg_RE(x));% ... %
%     @(x) myMISARE.con_RE(x));
optprob.options.TolFun = 1e-4;
optprob.options.Hessian = {'lbfgs'  5};
optprob.options.Diagnostics = 'on';
optprob.options.Display = 'iter-detailed';
optprob.options.PlotFcns = {@optimplotfval @optimplotconstrviolation @optimplotstepsize @optimplotfirstorderopt};
[wout,fval,exitflag,output] = fmincon(optprob);

total = 1; % Number of random starts

%% Run with multiple restarts
% tp = 'con_'; %[]
% tp = 'con_c8_'; %[]
% tp = 'con_pw0_'; %[]
delta_f = 1e-5;

max_try = 1;

tp = 'con_bound_';

rng(1980 + r)
scurr = rng;

for rr = 1:total
    rng(scurr)
    
    % Random starting pointw
    w0 = [];
    W0 = cell(size(A));
%     l = sqrt(3./O.V);
    for ff = M
        [u, s, v] = svd(randn(size(A{ff}')),'econ');
        W0{ff} = u*v';
%         W0{ff} = l(ff)*10*(rand(size(O.W{ff})) - .5);
    end
    scurr = rng;
    w0 = stackW(W0);
    myMISA = MISA(w0, M, S, X);
    
    myMISAobjective = @(x) myMISA.objective(x);
    
    optprob = getoptprob4(w0, myMISAobjective);
    [wout,fval,exitflag,output] = fmincon(optprob);
    Wout = unstackW(wout,myMISA.M,myMISA.C,myMISA.V);
    wbest = wout;
    
    mISI = MISI(Wout,A,myMISA.S);
    %save(['..' f 'output' f 'S_case' cas '_A_' casA f 'r' num2str(r,'%03d') ...
    %        f tp f tag '_' tagA '_' num2str(rr,'%03d') '.mat'], ...
    %        'scurr', 'w0', 'wout', 'mISI', 'fval', 'exitflag', 'output')
end

%% IVA Tulay
tp = 'tuln_noinit_';
rng(1980 + r)
scurr = rng;

for rr = 1:total
    rng(scurr)
    
    % Random starting pointw
    w0 = [];
    W0 = cell(size(A));
    % l = sqrt(3./V);
    for ff = M
        [u, s, v] = svd(randn(size(A{ff}')),'econ');
        W0{ff} = u*v';
        %     W0{ff} = l(ff)*2*(rand(size(A{ff}')) - .5);
    end
    scurr = rng;
    w0 = stackW(W0);
    
    in = cat(3,Xdat{M});
    W0iva = cat(3,W0{M});
    if strcmpi(tp, 'tuln_')
        error('should not go in here!!')
        [Wiva] = icatb_iva_second_order(in,'W_init',W0iva,'opt_approach','newton');
    elseif strcmpi(tp, 'tuln_noinit_')
        [Wiva] = icatb_iva_second_order(in,'opt_approach','newton');
    end
    %[Wiva] = icatb_iva_laplace(in,'maxIter',1024,'initW',Wiva,'whiten',false,'alpha0',1,'terminationCriterion','ChangeInCost','verbose',true);
    [Wiva] = icatb_iva_laplace(in,'maxIter',1024,'initW',Wiva,'whiten',false,'alpha0',1,'terminationCriterion','ChangeInW','verbose',true);
    for mm = 1:length(M)
        Wout{M(mm)} = Wiva(:,:,mm);
    end
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
% tp = 'tuln_';
% tp = 'tuln_noinit_';
R = 6;
for r = 1:R
    group = [group r*ones(1,total)];
    for rr = 1:total
        try
            load(['..' f 'output' f 'S_case' cas '_A_' casA f 'r' num2str(r,'%03d') ...
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

figure('position',[50   50   750   600])
hold on
B = 15;
owd = -50;
imhist_ = zeros(B,R);
for r = 1:R
%     [c, e] = histme_1D(mISIlist(group == r),[min(mISIlist)-eps max(mISIlist)+eps],B);
    ylim_ = min(owd, floor(min(mISIlist(:))/10)*10);
    [c, e] = histme_1D(mISIlist(group == r),[ylim_ 0],B);
    imhist_(:,r) = c;
end
es = e(1:B) + (e(2)-e(1))/2;
imagesc([1 R], [es(1) es(B)],imhist_,[0,total])
axis square tight
axis xy
ht = hot;
colormap(ht(64:-1:17,:))
set(gca,'xtick',1:R)
switch lower(cas)
    case '12'
        set(gca,'xticklabel',{0,.11,.23,.39,.5,.65})
          
end
set(gca,'ylim',[e(1) e(end)])
yticks = 0:-5:ylim_;
set(gca,'ytick',yticks(end:-1:1))
% set(gca,'ytick',e)
% set(gca,'yticklabel',cellstr(num2str(e,'%.1f')))
switch lower(tp)
    case 'con_bound_'
        title('MISA (with SOS) solving an IVA problem')
    case 'tuln_'
        title('IVA-GL solving an IVA problem')
    case 'tuln_noinit_'
        title('IVA-GL (self init.) solving an IVA problem')
end
ch = colorbar;
ylabel(ch, 'Counts')
switch lower(cas)
    case '12'
        xlabel('Maximum of correlation in subspaces')
          
end
ylabel('joint ISI [dB]')

hold on
for r = 1:R
    plot(group(group == r) + .6*rand(1,sum(group==r))-.3,mISIlist(group==r),'.k')
end

