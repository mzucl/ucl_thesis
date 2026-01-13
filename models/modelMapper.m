% RULE: Only the first function in the file (here: `modelMapper`) can be called from outside the file.
% All others (like `qZcUpdate`) are local to the file and accessible only from within it.

function modelMapper(data, ~, intermKVStore, params)
    Kc = params.K;
    Nc = size(data, 1);

    Zc = GaussianContainer("DS", Nc, true, Kc, zeros(Kc, 1));

    X = table2array(data);
    Nc = size(X, 1);
    Xc = ViewHandler(X, true);
    

    % === STEP 1: qZ update ===
    qZcUpdate(Zc, Xc, params);

    % === STEP 2: emit sufficient statistics ===
    stats.N = Nc;
    stats.E_ZZt = Zc.E_XXt;
    stats.X_times_E_Zt = Xc.X * Zc.E';
    stats.E_Z_times_Xt = Zc.E * Xc.X';
    stats.H_Z = Zc.H;
    stats.E_LnP_Z = Zc.E_LnP; 
    stats.X_col_sum = sum(Xc.X, 2);
    stats.Z_col_sum = sum(Zc.E, 2);
    stats.Tr_XtX = Xc.Tr_XtX;
        
    add(intermKVStore, 'stats', stats);
end

function qZcUpdate(Zc, Xc, params)
    W_    = params.W;
    mu_   = params.mu;
    tau_  = params.tau;
    it_   = params.it;
    covZ_ = params.covZ;

    tauExp = LogicUtils.ternary(it_ == 1, tau_.getExpInit(), tau_.E);
    muExp = LogicUtils.ternary(it_ == 1, mu_.getExpInit(), mu_.E);

    muNew = tauExp * covZ_ * W_.E_Xt * (Xc.X - muExp);
    Zc.updateDistributionsParameters(muNew, covZ_);
end