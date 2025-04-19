function data = GFA_Syn_2G()
    % Spherical noise - check the definition for tau!
    %
    %
    N_train = 400; 
    N_test = 100; % For predicting views
    N = N_train + N_test;

    M = 2;
    D = [50, 30]; % Dimensions of the views

    % Latent factors
    K = 4;
    Z = zeros(K, N);

    % Generate latent factors - row by row (in total K rows)
    n = 1:N;
    Z(1, :) = sin((n)/(N/20));
    Z(2, :) = cos((n)/(N/20));
    Z(3, :) = 2 * ((n)/N-0.5); 
    Z(4, :) = normrnd(0, 1, [N, 1]);

    % Noise precisions (tau) for each view; Spherical noise;
    tau = {5; 10};
    
    % Alpha (alpha) for each view K values
    alpha = zeros(M, K);
    alpha(1, :) = [1e6, 1, 1e6, 1];
    alpha(2, :) = [1, 1, 1, 1e6];
    
    % W matrix for each view
    W = cell(M, 1);

    % Train and test observarions: for all views M
    X_tr = cell(M, 1);
    X_te = cell(M, 1);
    
    for m = 1:M
        Dm = D(m);
        W{m} = zeros(Dm, K);
        for k = 1:K
            % Generate w_k (kth column of W) from p(W | alpha)
            alpha_m_k = alpha(m, k); % alpha (precision) for the view 'm' and column 'k'

            % normrnd(MU,SIGMA) returns an array of random numbers chosen from a
            % normal distribution with mean MU and standard deviation SIGMA.
            W{m}(:, k) = normrnd(0, 1/sqrt(alpha_m_k), [Dm, 1]);
        end

        % Generate X for view 'm'
        X = W{m} * Z + normrnd(0, 1/sqrt(tau{m}), [Dm, N]);

        % Get training and test data
        X_tr{m} = X(:, 1:N_train);
        X_te{m} = X(:, N_train+1:end);
    end
    
    % Latent variables for training the model    
    Z = Z(:, 1:N_train);
    
    % Store data and model parameters            
    data.X_tr = X_tr;
    data.X_te = X_te;
    data.W = W;
    data.Z = Z;
    data.tau = tau;
    data.alpha = alpha;
    data.K = K;
end
