function W = SNMTFFroInfer(V, opts),

Winit = opts.W0;
Hinit = opts.H0;

W = Winit;
H = Hinit;
k = size(Winit, 2);

W=W';

maxNumberThreads = opts.maxNumberThreads + 0.0;

W(:,:) = 0;

size(H)
size(V)

HH = H*H'; %(kxm)*(mxk)
VH = -H*double(V'); %(kxm)*(mxn)]
if isfield(opts, 'mu_2'),   VH = VH + opts.mu_2;  end
if isfield(opts, 'beta_2'),   HH = HH + 2 * opts.beta_2 * eye(k); end    
[W, ~] = doSNQPInfer(HH, VH, W, 1000000, 1e-9, 0, maxNumberThreads);

W = W';
end