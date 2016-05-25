addpath (genpath('..'));
clear
rand('seed',1);

%%% set the size/sparsity of the data
m = 10000;
n = 500;
k = 30;
rate = 0.3;

V = rand(n, m) * 100;

opts = getOptions();
opts.maxIter = 200;
opts.maxNumberThreads = 1;
maxNumCompThreads(opts.maxNumberThreads)
%% running GCD with trace=0, objGCD and timeGCD will just contain the final objective function and total running time
opts.W0 = rand(n, k) * max(max(V));
opts.H0 = rand(k, m);
[W, H, HIS] = SNMFFro(V, opts);
getSparsity(W)
getSparsity(H)

%norm(V -W*H, 'fro')^2/2.0;
opts.W0 = rand(n, k) * max(max(V));
opts.H0 = rand(k, m);
[W, H, HIS] = NMFFro(V, opts);
getSparsity(W)
getSparsity(H)


