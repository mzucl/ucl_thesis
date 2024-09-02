function data = generateTwoViews()
    Ntrain = 400; Ntest = 100;
    N = Ntrain + Ntest;
    M = 2;
    D = [50, 30]; % Dimensions of the views

    % Latent factors
    K = 4;
    Z = zeros(K, N);

    x = 1:N;

    Z(1, :) = sin(x / (N/20));
    Z(2, :) = cos(x / (N/20));
    Z(3, :) = 2 * (x/N - 0.5); 
    Z(4, :) = normrnd(0, 1, [N, 1]);          
    
    % Specify noise precisions manually
    tau = cell(1, length(D));
    tau{1} = 5 * ones(1, D(1));
    tau{2} = 10 * ones(1, D(2));
    
    % Specify alphas manually
    alpha = zeros(M, K);
    alpha(1, :) = [1e6, 1, 1e6, 1];
    alpha(2, :) = [1,   1, 1,   1e6];
    
    W = cell(1, M);
    X_train = cell(1, M);
    X_test = cell(1, M);
    
    for i = 1:length(D)
        W{i} = zeros(D(i), K);
        for t = 1:K
            % generate W from p(W|alpha)
            W{i}(:,t) = normrnd(0, 1./sqrt(alpha(i,t)), [D(i), 1]);
        end
        X = zeros(N, D(i));
        for j = 1:D(i)
            % generate X from the generative model
            X(:,j) = Z * W{i}(j,:)' + normrnd(0, 1./sqrt(tau{i}(j)), [N, 1]);
        end
        % Get training and test data
        X_train{i} = X(1:Ntrain,:); % Training data
        X_test{i} = X(Ntrain+1:end,:); % Test data
    end
    
    % latent variables for training the model    
    Z = Z(1:Ntrain,:);
    
    % Store data and model parameters            
    data.X_tr = X_train;
    data.X_te = X_test;
    data.W = W;
    data.Z = Z;
    data.tau = tau;
    data.alpha = alpha;
    data.trueK = K;
end
