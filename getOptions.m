function opts = getOptions(),
    opts.type = 'PLAIN';
    opts.maxIter = 10;
    opts.tolerance = 1e-3;
    opts.threads = 1;
    opts.verbose = 1;
    opts.maxNumberThreads = 1;
    opts.params = [0, 0, 0, 0]';
    opts.beta_2 = 0.0;
end