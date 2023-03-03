function myMISAKRE = setup_basic_MISAKRE(sim1, beta, gradtype, sc, RElambda)

ut = utils;

preX = false;

beta = beta*ones(size(sim1.S{sim1.M(1)},1),1);
eta = ones(size(sim1.S{sim1.M(1)},1),1);
lambda = ones(size(sim1.S{sim1.M(1)},1),1);
% gradtype = 'regular';

REtype = 'NMSE';
REapproach = 'PINV';
% RElambda = 1e-2;
REref = {}; %REref = sim1.ref;
REreftype = 'linearreg';
REreflambda = {.9};
rC = sim1.predictors;

myMISAKRE = MISAKRE(ut.stackW(sim1.Wi), ...
    sim1.M, sim1.S, sim1.genX(), beta, eta, lambda, gradtype, sc, preX, ...
    REtype, REapproach, RElambda, ...
    REref, REreftype, REreflambda, rC);