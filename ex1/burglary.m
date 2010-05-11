

N   = 4;
dag = zeros(N,N);
B = 1; E = 2; A = 3; R = 4;
dag(B, A) = 1;
dag(E, [A R]) = 1;

discrete_nodes = 1:N;
node_sizes = 2*ones(1,N);

bnet = mk_bnet(dag, node_sizes, 'observed', [R A]);
bnet.CPD{B} = tabular_CPD(bnet, B, [0.99 0.01]);
bnet.CPD{E} = tabular_CPD(bnet, E, [1-10^-6 10^-6]);
bnet.CPD{R} = tabular_CPD(bnet, R, [1 0 0 1]);
bnet.CPD{A} = tabular_CPD(bnet, A, [1-0.001 1-0.95 1-0.41 1-0.98 0.001 0.95 0.41 0.98]);

engine = jtree_inf_engine(bnet);

evidence = cell(1,N);
evidence{A} = 2;
evidence{R} = 2;

engine
[engine, loglik] = enter_evidence(engine, evidence);
engine

marg = marginal_nodes(engine, E);
marg.T
