% RULE: Only the first function in the file (here: `modelMapper`) can be called from outside the file.
% All others (like `qZcUpdate`) are local to the file and accessible only from within it.

function modelMapper(data, info, intermKVStore, params)
    [~, name, ~] = fileparts(info.Filename);   % Get filename without path/extension
    tokens = regexp(name, 'Xdata_part(\d+)', 'tokens');
    chunkId = str2double(tokens{1}{1});
    

    % ZcPath = fullfile("ZcStorage", sprintf("Zc_%d.mat", chunkId));

    W_ = params.W;
    mu_ = params.mu;

    Kc = params.K;
    Nc = size(data, 1);
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

    Zc = GaussianContainer("DS", Nc, true, Kc, zeros(Kc, 1));

    X = table2array(data);
    Nc = size(X, 1);
    Xc = ViewHandler(X, true);
    

% If you prefer calling save(ZcPath, 'Zc') directly, you could subclass from matlab.mixin.CustomSaveLoad to integrate with MATLAB's serialization system more natively (let me know if you want that).
% 
% Saving private properties is optional—do it only if needed for correctness or performance.
% 



    % if isfile(ZcPath)
    %     % ZcStruct = load(ZcPath);
    %     % Zc = ZcStruct.Zc;
    % 
    % 
    %     load(ZcPath, 'ZcStruct');
    %     Zc = GaussianContainer.loadobj(ZcStruct);  % reconstruct object
    % else
    %     initZMu = randn(Kc, 1);
    %     Zc = GaussianContainer("DS", Nc, true, Kc, initZMu);
    % end


    % === STEP 1: qZ update ===
    qZcUpdate(Zc, Xc, params);

    % === Save updated Zc back to file ===
    % if ~isfolder("ZcStorage")
    %     mkdir("ZcStorage");
    % end
    % save(ZcPath, 'Zc');
    % 
    % ZcStruct = Zc.saveobj;       % convert handle object to struct
    % save(ZcPath, 'ZcStruct');     % save struct to file

    % === STEP 2: emit sufficient statistics ===
    stats.N = Nc;
    stats.E_ZZt = Zc.E_XXt;




    stats.X_times_E_Zt = Xc.X * Zc.E';

    



    % FINE!!!
    stats.E_Z_times_Xt = Zc.E * Xc.X';
    stats.E_Z = Zc.E;
    stats.H_Z = Zc.H;
    stats.E_LnP_Z = Zc.E_LnP; 
    stats.X_col_sum = sum(Xc.X, 2);
    stats.Tr_XtX = Xc.Tr_XtX;
        
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