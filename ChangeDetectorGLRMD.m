classdef ChangeDetectorGLRMD < ChangeDetectorGLR
    % Multi-dimensional implementation of GLR
    
    properties(Access = protected)
        % The number of feature extraction methods used. Must be >1.     
        featureExtractorsCount;
    end % class data members
    
    methods(Access = public)
        
        function obj = ChangeDetectorGLRMD()
            % Constructor
            obj@ChangeDetectorGLR();
        end % class constructor

        function [obj, allScores, windowsTime, features] = detectChange(obj, signal) 
            % 'signal' is the signal we extract the features from.
            % It can be 3D: [signalsCount, electrodesCount, samplesCount].
            % It can be 2D: [electrodesCount, samplesCount].
            % 3D format is recommended though! 
            % The number of samples must be greater or equal to 2 but 
            % generally we need far more that. 
            % If there is more than one target, 'allScores' will be 3D, 
            % where the first dimension indexes the signal, the second
            % dimension indexes the target, and the third dimension indexes
            % the score. This happens because we do not know the correct
            % target. Thus, the algorithm assumes each target is the
            % correct target and performs the computation bases on that
            % assumption. Thus for N target, we get N different results for
            % each signal.   
            
            [obj, signal] = setOtherParameters(obj, signal);            
            obj.featureExtractorsCount = length(obj.featureExtractor);
                        
            % Creat a copy of all feature extraction method lest we lose
            % them. We process each method one at a time. 
            allFeatureExtractors = obj.featureExtractor;

            isMultiTargets = zeros(obj.featureExtractorsCount, 1);
            
            % features will have size (feature dim, batch, features)
            featuresSingleTarget = [];
            
            % This will have size (targets, batch, features)
            featuresMultiTarget = [];
                       
            for i = 1:obj.featureExtractorsCount
                
                obj.featureExtractor = allFeatureExtractors{i};
                
                % Move the signal from time domain to feature space
                [obj, windowsTime] =...
                    extractFeatureSignal(obj, signal);
                
                % Hopefully, features must be all real but for short window
                % lengths some feature extraction methods may fail.
                % I should look into this sometime.
                obj.featureSignal = real(obj.featureSignal);
                obj.featureSignal = normalize(obj.featureSignal');
                obj.featureSignal = obj.featureSignal';
                
                if isa(allFeatureExtractors{i}, ...
                        'FeatureExtractionTemplateMatching') &&...
                        getTargetsCount(allFeatureExtractors{i}) > 1                
                    
                    isMultiTargets(i) = 1;
                    targetsCount = ...
                        getTargetsCount(allFeatureExtractors{i});  
                    
                    obj.featureSignal = reshape(obj.featureSignal,...
                        [1, obj.signalsCount, targetsCount,...
                        size(obj.featureSignal, 2)]);
                   
                    obj.featureSignal = permute(obj.featureSignal,...
                        [1, 3, 2, 4]);
                    
                    if isempty(featuresMultiTarget)
                        featuresMultiTarget = obj.featureSignal;
                    else
                        featuresMultiTarget = cat(1, ...
                            featuresMultiTarget, obj.featureSignal);
                    end
                    
                    continue;                    
                end

                obj.featureSignal = reshape(obj.featureSignal,...
                    [obj.signalsCount,...
                    size(obj.featureSignal, 1)/obj.signalsCount,...
                    size(obj.featureSignal, 2)]);
                
                obj.featureSignal = permute(obj.featureSignal, [2, 1, 3]);
                
                if isempty(featuresSingleTarget)
                    featuresSingleTarget = obj.featureSignal;
                else
                    featuresSingleTarget = cat(1,...
                        featuresSingleTarget, obj.featureSignal);
                end             
                
            end % for i
                        
            if ~isempty(featuresMultiTarget)
                
                % Put the target dimension first, so we can squeeze it.
                featuresMultiTarget = ...
                    permute(featuresMultiTarget, [2, 1, 3, 4]);
                
                allScores = zeros(targetsCount, obj.signalsCount,...
                    size(obj.featureSignal, 2));
            
                for target = 1:targetsCount
                    features = squeeze(featuresMultiTarget(1, :, :, :));
                    
                    if ismatrix(features)
                        features = reshape(features, [1, size(features)]);
                    end
                    
                    features = cat(1, features, featuresSingleTarget);
                    obj.featureSignal = permute(features, [2, 1 , 3]);
                    scores = computeAllPossibleChangeTimes(obj);   
                    
                    if target == 1
                        allScores = zeros([targetsCount, size(scores)]);
                    end 
                    
                    allScores(i, :) = scores(:);
                end % for target
                
                % Batch signals first, then targets, then likelihoods
                allScores = permute(allScores, [2, 1, 3]);
                
            else
                obj.featureSignal = permute(featuresSingleTarget, [2, 1 , 3]);
                allScores = computeAllPossibleChangeTimes(obj); 
            end                
            
%             plotLikelihoods(obj, allScores(1:10:50, :), squeeze(obj.featureSignal(1, 1, :)), windowsTime);

              features = obj.featureSignal;
        end % function detectChange
        
    end % class public services
    
    methods(Access = protected)
                
        function likelihoodSum = ...
                computeCumulativeLikelihood(obj, j, k)
            % Compute cumulative some of likelihood ratios for all samples
            % between j and k. Assume normal distribution. Assume
            % independence of samples. 
                       
            if k > size(obj.featureSignal, 3)
                error('Stop index must not exceed features count.');
            end
            
            if j > k
                error('Window length is negative (wrong indexing)');
            end
            
            if j == 1
                estimatedMean0 = obj.featureSignal(:, :, 1);
                estimatedMean1 = mean(obj.featureSignal(:, :, 2:k), 3);
            else
                estimatedMean0 = mean(obj.featureSignal(:, :, 1:j-1), 3);
                estimatedMean1 = mean(obj.featureSignal(:, :, j:k), 3);
            end
            
            covarianceMatrix = zeros(size(obj.featureSignal, 1), ...
                size(obj.featureSignal, 2), ...
                size(obj.featureSignal, 2));
            
            % Use the larger portion of the signal to compute the variance.
            % Because we assume equal variances, the numerical value of the
            % variance does not matter in selecting the hypothesis.
            % Nonetheless, it helps us get more reallistic likelihoods. 
            for i = 1:size(obj.featureSignal, 1)
                
                if j-1 > k/2
                    covarianceMatrix(i, :, :) = ...
                        cov(squeeze(obj.featureSignal(i, :, 1:j-1))');
                else
                    covarianceMatrix(i, :, :) = ...
                        cov(squeeze(obj.featureSignal(i, :, j:k))');
                end 
                                
            end % for i
            
            middleSample = 0.5 * (k - j + 1);

            likelihoodSum = computeLikelihoodsNormal(obj,...
                estimatedMean0, estimatedMean1, covarianceMatrix, ...
                middleSample);
                      
        end % function computeCumulativeLikelihood 
        
        function likelihood = ...
                computeLikelihoodsNormal(obj,...
                mean0, mean1, covarMatrix, middleSample)
            % Compute the likelihood ratio for sample yi
            % Assume yi has normal distributions for both hypotheses
            % Assume both hypotheses have the same variance.
            
            deltaMean = (mean1 - mean0);
            likelihood = zeros(size(mean1, 1), 1);
            
            for i = 1:size(mean0, 1)
                meanDifference = deltaMean(i, :);
                currentCovariance = squeeze(covarMatrix(i, :, :));
                loglike = middleSample * meanDifference *...
                    inv(currentCovariance) * meanDifference';
                likelihood(i) = loglike;
            end % for i
            
        end % end function compute likelihoodNormal
        
        function plotLikelihoods(obj, s, fs, wt)
            figure;
            plot(wt(3, :), s', 'LineWidth', 3);
            xlabel('Time (s)');
            ylabel('Log likelihood');
            grid on;
            myLegend = legend('T1', 'T11', 'T21', 'T31', 'T41');
            myLegend.NumColumns = 2;
            set(gca, 'FontSize', 40);
            ylim([1.3*min(s(:)), 1.3*max(s(:))]);
            xlim([0, max(wt(3, :))]);
            
            figure;
            plot(wt(2, :), squeeze(obj.featureSignal(3, 1, :)), 'LineWidth', 3);
            xlabel('Time (s)');
            ylabel('Feature Amplitude');
            grid on;
%             myLegend = legend('E14', 'E15', 'E16');
%             myLegend.NumColumns = 2;
            set(gca, 'FontSize', 40);
            ylim([1.3*min(fs(:)), 1.3*max(fs(:))]);
            xlim([0, max(wt(3, :))]);
            
        end % function plotLikelihoods
        
    end % class utilities
    
end % class ChangeDetectorGLRMD