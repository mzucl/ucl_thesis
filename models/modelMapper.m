% RULE: Only the first function in the file (here: `modelMapper`) can be called from outside the file.
% All others (like `qZcUpdate`) are local to the file and accessible only from within it.

function modelMapper(data, ~, intermKVStore, params)
    persistent Zc;
    persistent Xc;

    W_ = params.W;
    mu_ = params.mu;

    Kc = params.K;
    Nc = size(data, 1);

    if isempty(Zc)
        initZMu = randn(Kc, 1);
        %                      type, size_, cols, dim,   mu
        Zc = GaussianContainer("DS",  Nc,   true,  Kc, initZMu);
    end

    if isempty(Xc)
        Xc = ViewHandler(table2array(data), true);
    end


    % STEP 1: qZUpdate
    % TODO: Covariance is shared, don't recompute it!
    qZcUpdate(Zc, Xc, params);

    % STEP 2: Compute sufficient stats for reducer
    stats.N = Nc;
    stats.E_ZcZct = Zc.E_XXt;
    stats.Tr_XctXc = Xc.Tr_XtX;
    stats.Xc_col_sum = sum(Xc.X, 2);
    stats.E_WZc_col_sum = sum(W_.E * Zc.E, 2);
    stats.E_Z_times_centered_data_t = Zc.E * (Xc.X' - mu_.E_Xt);

    expWZ = W_.E * Zc.E;
    stats.Tr_E_WZ_time_Xt = sum(Xc.X(:) .* expWZ(:));

    stats.H = Zc.H;
    stats.E_LnP = Zc.E_LnP; 
        
    add(intermKVStore, 'stats', stats);
end

function qZcUpdate(Zc, Xc, params)
    W_ = params.W;
    mu_ = params.mu;
    tau_ = params.tau;
    it_ = params.it;
    K_ = params.K;

    tauExp = Utility.ternary(it_ == 1, tau_.getExpInit(), tau_.E);
    muExp = Utility.ternary(it_ == 1, mu_.getExpInit(), mu_.E);

    covNew = Utility.matrixInverse(eye(K_) + tauExp * W_.E_XtX);
    muNew = tauExp * covNew * W_.E_Xt * (Xc.X - muExp);

    Zc.updateDistributionsParameters(muNew, covNew);
end