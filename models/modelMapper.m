% RULE: Only the first function in the file (here: `modelMapper`) can be called from outside the file.
% All others (like `qZcUpdate`) are local to the file and accessible only from within it.

function modelMapper(data, info, intermKVStore, params)
    [~, name, ~] = fileparts(info.Filename);   % Get filename without path/extension
    tokens = regexp(name, 'Xdata_part(\d+)', 'tokens');
    chunkId = str2double(tokens{1}{1});
    

    ZcPath = fullfile("ZcStorage", sprintf("Zc_%d.mat", chunkId));

    W_ = params.W;
    mu_ = params.mu;

    Kc = params.K;
    % Nc = size(data, 1);
    % 
    % if isempty(Zc)
    %     initZMu = randn(Kc, 1); % I think  I don;t need this init, becuase it is update below right away
    %     %                      type, size_, cols, dim,   mu
    %     Zc = GaussianContainer("DS",  Nc,   true,  Kc, initZMu);
    % end
    % 
    % if isempty(Xc)
    %     Xc = ViewHandler(table2array(data), true);
    % end


    X = table2array(data);
    Nc = size(X, 1);
    Xc = ViewHandler(X, true);
    
    if isfile(ZcPath)
        ZcStruct = load(ZcPath);
        Zc = ZcStruct.Zc;
    else
        initZMu = randn(Kc, 1);
        Zc = GaussianContainer("DS", Nc, true, Kc, initZMu);
    end


    % === STEP 1: qZ update ===
    qZcUpdate(Zc, Xc, params);

    % === Save updated Zc back to file ===
    if ~isfolder("ZcStorage")
        mkdir("ZcStorage");
    end
    save(ZcPath, 'Zc');

    % === STEP 2: emit sufficient statistics ===
    stats.N = Nc;
    stats.E_ZcZct = Zc.E_XXt;
    stats.E_Zc_times_centered_data_t = Zc.E * (Xc.X' - mu_.E_Xt);

    stats.Xc_col_sum = sum(Xc.X, 2);
    stats.E_WZc_col_sum = sum(W_.E * Zc.E, 2);

    stats.Tr_XctXc = Xc.Tr_XtX;
   
    expWZ = W_.E * Zc.E;
    stats.Tr_E_WZ_times_Xt = sum(Xc.X(:) .* expWZ(:));

    stats.H = Zc.H;
    stats.E_LnP = Zc.E_LnP; 
        
    add(intermKVStore, 'stats', stats);
end

function qZcUpdate(Zc, Xc, params)
    W_ = params.W;
    mu_ = params.mu;
    tau_ = params.tau;
    it_ = params.it;
    covZ_ = params.covZ;

    tauExp = Utility.ternary(it_ == 1, tau_.getExpInit(), tau_.E);
    muExp = Utility.ternary(it_ == 1, mu_.getExpInit(), mu_.E);
    muNew = tauExp * covZ_ * W_.E_Xt * (Xc.X - muExp);

    Zc.updateDistributionsParameters(muNew, covZ_);
end