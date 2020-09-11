classdef ChangeDetectorGLR < ChangeDetector
    % A class for detecting changes using generalized likelihood ratio
    % algorithm. The class assumes that the samples are iid and have a
    % normal distribution with the same variance before and after change
    % but different mean values. This class only works with one-dimensional
    % signals (in the feature domain. The signal can be still
    % multi-dimensional if time domain. For example, CCA works on 
    % multi-dimensional signals in time domain but the generated features 
    % are scalers). 
    
    properties(Access = protected)
        % Empty
    end % class protected data members
    
    methods(Access = public)
        
        function obj = ChangeDetectorGLR()
            % Constructor
            
            obj@ChangeDetector();            
        end % class constructor
        
        function [obj, allScores, windowsTime, features] = detectChange(obj, signal)
            % 'signal' is the signal we extract the features from.
            % It can be 3D: [signalsCount, electrodesCount, samplesCount].
            % It can be 2D: [electrodesCount, samplesCount].
            % The number of samples must be greater or equal to 2.              
            % If there is more than one target, 'allScores' will be 3D, 
            % where the first dimension indexes the signal, the second
            % dimension indexes the target, and the third dimension indexes
            % the score. This happens because we do not know the correct
            % target. Thus, the algorithm assumes each target is the
            % correct target and performs the computation bases on that
            % assumption. Thus for N target, we get N different results for
            % each signal. 
              
            obj = setOtherParameters(obj, signal);
            
            % Move the signal from time domain to feature space
            [obj, windowsTime] =...
                extractFeatureSignal(obj, signal);
            
            % Test all hypotheses and compute the log likelihood for each
            % of them.
            allScores = computeAllPossibleChangeTimes(obj);       
            
            % If multi-target, reshape the data for easier access
            if obj.featuresMatrixDim(2) > 1
                targetsCount = size(allScores, 2);
                allScores = reshape(allScores, [obj.signalsCount,...
                    obj.featuresMatrixDim(2), targetsCount]);
                obj.featureSignal = reshape(obj.featureSignal,...
                    [obj.signalsCount, obj.featuresMatrixDim(2), ...
                    targetsCount]);
            end 
            
            features = obj.featureSignal;
%             plotLikelihoods(obj, allScores(1:10:50, :), obj.featureSignal(1:10:50, :), windowsTime);
        end % function detectChange

    end % class public services
    
    methods(Access = protected)

        function allScores = ...
                computeAllPossibleChangeTimes(obj)
            % For a process with K samples, we test the following 
            % hypotheses. First, we assume that the change happens at the
            % first sample (in theory we never accept this hypothesis. I
            % have added it for simpler programming.) Then we test the
            % hypothesis of change happening at the second sample. Then,
            % the third sample, and so on.
            
            % Extract the size of last dimension, which holds the samples.
            k = size(obj.featureSignal);
            k = k(end);
            
            % Save the log likelihood of each hypothesis.
            % We generally pick the hypothesis with the largest log
            % likihood. We assume each sample can be the moment where the
            % change happens. Thus, we need to test k hypotheses. 
            allScores = zeros(size(obj.featureSignal, 1), k);
            
            if k <= 1
                error('process too short.\n');
            end
            
            for j = 1:k
                score = computeCumulativeLikelihood(obj, j, k);
                allScores(:, j) = score;
            end % for j
            
        end % end function computeAllPossibleChangeTimes
        
        function likelihoodSum = ...
                computeCumulativeLikelihood(obj, j, k)
            % Compute cumulative some of likelihood ratios for all samples
            % between j and k. Assume normal distribution. Assume
            % independence of samples. 
            
            % Ensure we do not exceed the array.
            if k > size(obj.featureSignal, 2)
                error('stop index must not exceed features count.\n');
            end
            
            % This condition should never happen.
            % If it ever does, there must be an error in computing j and k.
            if j > k
                error('Window length is negative (wrong indexing)');
            end
            
            % Typically, we assume that change cannot happen at the very
            % first sample. However, I created a condition to process this
            % situation gracefully. The computed results for when we assume
            % that the change happens at the very first sample should
            % be discarded because the theory does not account for it. 
            if j == 1
                estimatedMean0 = obj.featureSignal(:, 1);
                estimatedMean1 = mean(obj.featureSignal(:, 2:k), 2);
            else
                estimatedMean0 = mean(obj.featureSignal(:, 1:j-1), 2);
                estimatedMean1 = mean(obj.featureSignal(:, j:k), 2);
            end
            
            % Use the larger portion of the signal to compute the variance.
            % Because we assume equal variances, the numerical value of the
            % variance does not matter in selecting the hypothesis.
            % Nonetheless, it helps us get reallistic likelihoods. 
            if j-1 > k/2
                variance = var(obj.featureSignal(:, 1:j-1), 0, 2);
            else
                variance = var(obj.featureSignal(:, j:k), 0, 2);
            end 
            
            % Compute the likelihood ratios for each sample between sample
            % j and sample k. Then sum them to get the total. We add these
            % values because we assume independence. 
            allLikelihoods = computeLikelihoodsNormal(obj,...
                obj.featureSignal(:, j:k), ...
                estimatedMean0, estimatedMean1, variance);
            
            likelihoodSum = sum(allLikelihoods, 2);            
        end % function computeCumulativeLikelihood
        
        function likelihood = ...
                computeLikelihoodsNormal(obj, yi, mean0, mean1, variance)
            % Compute the likelihood ratio for sample yi.
            % Assume yi has normal distributions for both hypotheses.
            % Assume both hypotheses have the same variance.
            
            likelihood = ((yi - mean0).^2 - (yi - mean1).^2) ....
                ./ (2 * variance);
            
            likelihood = likelihood / 2;
        end % end function compute likelihoodNormal
        
        function plotLikelihoods(obj, s, fs, wt)
            figure;
            plot(wt(2, :), s', 'LineWidth', 3);
            xlabel('Time (s)');
            ylabel('Log likelihood');
            grid on;
            myLegend = legend('T1', 'T11', 'T21', 'T31', 'T41');
            myLegend.NumColumns = 2;
            set(gca, 'FontSize', 28);
            xlim([0, 7]);
            
            figure;
            plot(wt(2, :), fs, 'LineWidth', 3);
            xlabel('Time (s)');
            ylabel('Feature Amplitude');
            grid on;
            myLegend = legend('T1', 'T11', 'T21', 'T31', 'T41');
            myLegend.NumColumns = 2;
            set(gca, 'FontSize', 28);
            xlim([0, 7]);
        end % function plotLikelihoods
        
    end % class utilities   
    
end % class ChangeDetectorGLR