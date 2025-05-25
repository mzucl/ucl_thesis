function modelMapper(data, ~, intermKVStore, params)
    persistent Zc;
    persistent Xc;

    W_ = params.W;
    mu_ = params.mu;
    tau_ = params.tau;
    it_ = params.iter;

    Kc = W_.dim; % Same as obj.K.Val
    Nc = size(data, 2);

    if isempty(Zc)
        initZMu = randn(Kc, 1);
        %                      type, size_, cols, dim,   mu
        Zc = GaussianContainer("DS",  Nc,   true,  Kc, initZMu);
    end

    if isempty(Xc)
        Xc = ViewHandler(table2array(data), false);
    end


    % STEP 1: qZUpdate
    % TODO: Covariance is shared, don't recompute it! HOW???
    tauExp = Utility.ternary(it_ == 1, tau_.getExpInit(), tau_.E);
    muExp = Utility.ternary(it_ == 1, mu_.getExpInit(), mu_.E);

    covNew = Utility.matrixInverse(eye(Kc) + tauExp * W_.E_XtX);
    muNew = tauExp * covNew * W_.E_Xt * (Xc.X - muExp);

    Zc.updateDistributionsParameters(muNew, covNew);

    % STEP 2: Compute sufficient stats for reducer
    stats.E_ZcZct = Zc.E_XXt;
    stats.Tr_XctXc = Xc.Tr_XtX;
    stats.Xc_col_sum = sum(Xc.X, 2);
    stats.E_WZc_col_sum = sum(W_.E * Zc.E, 2);
    
    add(intermKVStore, 'stats', stats);
end