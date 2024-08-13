% TODO (medium): Create a base class for all models with optimization params
% and some other stuff
classdef GFA < handle
    properties
        K               % Number of latent dimensions/principal components
    
        Z               % [K x N] GaussianDistributionContainer [size: N; for each latent variable zn]
        Groups          % An array of GFAGroup instances

        % Optimization parameters
        maxIter
        tol
        % They all share Z, also it shouldn't be a copy it should be a
        % reference!!!
    end

    
    properties (Dependent)
        M       % Number of groups
    end
    
    methods
        %% Deep copy and operators overloading
        function newObj = copy(obj)
            newObj = GammaDistribution();
            
            newObj.a = obj.a;
            newObj.b = obj.b;
            
            % Copy the prior (manually)
            if ~Utility.isNaN(obj.prior)
                newObj.prior = GammaDistribution();
                newObj.prior.a = obj.prior.a;
                newObj.prior.b = obj.prior.b;
            end
        end

        function isEqual = eq(obj1, obj2)
            if ~Utility.isNaNOrInstanceOf(obj1, 'GammaDistribution') || ~Utility.isNaNOrInstanceOf(obj2, 'GammaDistribution')
                isEqual = false;
                return;
            end

            % Both are NaN
            if Utility.isNaN(obj1) && Utility.isNaN(obj2)
                isEqual = true;
                return;
            % Only one is NaN
            elseif xor(Utility.isNaN(obj1), Utility.isNaN(obj2))
                isEqual = false;
                return;
            end
                
            % Both are set - compare them!    
            isEqual = obj1.a == obj2.a && obj1.b == obj2.b;

            % If parameters a and b are not equal, they are different objects
            % regardless of the prior!
            %   Second part of the condition is added because obj1.prior ==
            %   obj2.prior won't call 'isEqual' if both are NaN, and we are
            %   relaying on that call for comparing priors.
            if ~isEqual || Utility.isNaN(obj1.prior) && Utility.isNaN(obj2.prior)
                return;
            end

            isEqual = obj1.prior == obj2.prior;
        end

        function isNotEqual = ne(obj1, obj2)
            isNotEqual = ~eq(obj1, obj2);
        end



        % [NOTE] We need to deal with the form of the datasets here (e.g.
        % is it a table (for now it is), or it is file path so we need to
        % import the data. Also, we should check if all datasets have the
        % same number of observations.
        %% Constructors
        function obj = GFA(X)
            
            switch nargin
                case 0
                    obj.a = Constants.DEFAULT_GAMMA_A;
                    obj.b = Constants.DEFAULT_GAMMA_B;

                case 1
                    % If the ONLY parameter is of a type GammaDistribution that is
                    % the prior and we initialize the obj and its prior
                    % using that value
                    if Utility.areAllInstancesOf(a, 'GammaDistribution')
                        obj = a.copy();
                        obj.prior = a.copy();

                    % The one parameter passed in is the value for 'a', set
                    % both 'a' and 'b' to that value
                    elseif Utility.isSingleNumber(a)
                        obj.a = a;
                        obj.b = a;
                    else
                        error(['##### ERROR IN THE CLASS ' class(obj) ': Invalid arguments passed.']);
                    end
            
                case {2, 3}
                    obj.a = a;
                    obj.b = b;
                    if nargin == 3 && ~Utility.isNaN(prior) % 'prior' is passed in
                        obj.prior = prior.copy();
                    end
                otherwise
                    error(['##### ERROR IN THE CLASS ' class(obj) ': Too many arguments passed into the constructor.']);
            end
            % Set initial expectation to the real expectation
            obj.setExpInit(obj.Expectation);
        end



        %% Update methods



        %% Getters
        function value = get.D(obj)
            value = gamrnd(obj.view.D);
        end
    end
end