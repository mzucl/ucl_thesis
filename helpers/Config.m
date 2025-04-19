classdef Config < matlab.mixin.SetGet
    % Define nested classes (these will represent groups of settings)
    classdef Distributions < handle
        properties
            Gaussian
            Gamma
        end

        methods
            function obj = Distributions()
                obj.Gaussian.mu = 0;
                obj.Gaussian.precision = 1;
                obj.Gaussian.dim = 1;

                obj.Gamma.a = 1e-14;
                obj.Gamma.b = 1e-14;
            end
        end
    end

    classdef Optimization < handle
        properties
            MaxIter
            Tol
        end

        methods
            function obj = Optimization()
                obj.MaxIter = 5000;
                obj.Tol = 1e-6;
            end
        end
    end

    classdef Experiment < handle
        properties
            stabilityRun
            modelSelectionIter
            latentFactorThreshold
        end

        methods
            function obj = Experiment()
                obj.stabilityRun = 2;
                obj.modelSelectionIter = 2;
                obj.latentFactorThreshold = 1e-6;
            end
        end
    end

    % Main class properties
    properties (Access = private)
        Distributions
        Optimization
        Experiment
        Debug
    end

    methods (Access = private)
        function obj = Config()
            % Initialize nested classes
            obj.Distributions = obj.Distributions();
            obj.Optimization = obj.Optimization();
            obj.Experiment = obj.Experiment();

            % Additional settings (Debug, etc.)
            obj.Debug.validate = true;
            obj.Debug.verbose = true;
            obj.Debug.epsilon = 1e-2;
        end
    end

    methods (Static)
        function obj = getInstance()
            persistent instance
            if isempty(instance)
                instance = Config();
            end
            obj = instance;
        end
    end
end
