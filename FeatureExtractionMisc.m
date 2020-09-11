classdef FeatureExtractionMisc < FeatureExtractionMethod
    
    properties(Access = protected)
        featureType;
        featureIsMean;
        featureIsVariance;
        featureIsEnergy;
        featureIsKurtosis;
        featureIsRegression;
    end % class datamembers 
    
    methods(Access = public)
        
        function obj = FeatureExtractionMisc()
            % Class constructor
            obj@FeatureExtractionMethod();
            obj.featureType = 0;
            obj.featureIsMean = 0;
            obj.featureIsVariance = 1;
            obj.featureIsEnergy = 2;  
            obj.featureIsKurtosis = 3;
            obj.featureIsRegression = 4;
        end % function FeatureExtractionMisc
                
    end % class public services
    
    methods(Access = protected)
        
        function testvalues = computeFeatures(obj, Y)
            
            if obj.featureType == obj.featureIsMean
                testvalues = computeFeaturesMean(obj, Y);
                return;
            end 
            
            if obj.featureType == obj.featureIsVariance
                testvalues = computeFeaturesVariance(obj, Y);
                return;
            end
                      
            if obj.featureType == obj.featureIsEnergy
                testvalues = computeFeaturesEnergy(obj, Y);
                return;
            end
            
            if obj.featureType == obj.featureIsKurtosis
                testvalues = computeFeaturesKurtosis(obj, Y);
                return;
            end 
            
            if obj.featureType == obj.featureIsRegression
                testvalues = computeFeaturesRegression(obj, Y);
                return;
            end 
            
            error('Undefined feature extraction method.');
            
        end % class computeFeatures
        
        function testvalues = computeFeaturesMean(obj, Y)
            % Return the mean value of each signal
            
            testvalues = transpose(mean(Y, 1));
            
        end % function computeFeaturesMean
        
        function testvalues = computeFeaturesVariance(obj, Y)
            % Return the variance of each signal 
            
            testvalues = transpose(var(Y, 1));
            
        end % function computeFeaturesVariance
        
        function testvalues = computeFeaturesEnergy(obj, Y)
            % Return the energy of the signal
            
            testvalues = transpose(sum(Y.^2, 1));
            
        end % function computeFeaturesEnergy
        
        function testvalues = computeFeaturesKurtosis(obj, Y)
            
            correctBiasFlag = 0;
            testvalues = transpose(kurtosis(Y, correctBiasFlag, 1));
            
        end % function computeFeaturesKurtosis
        
        function testvalues = computeFeaturesRegression(obj, signal)
            
            targetFrequency = 10;
            testvalues =  regressionSineKernel(signal,...
                obj.samplingFrequency, targetFrequency);
            
        end % function computeFeaturesRegression
        
        function obj = performInitializations(obj)
            % Empty
            
        end % function performOtherPrecomputations     
        
        function obj = parseInputVariables(obj, options)

            obj = parseInputVariables@FeatureExtractionMethod(...
                obj, options); 
            
            parser = inputParser;
            parser.KeepUnmatched = true;
            addParameter(parser, 'featureType', 'mean');
            parse(parser, options{:});

            obj = setFeatureType(obj, parser.Results.featureType);

        end % function parseInputVariables
        
        function obj = setFeatureType(obj, type)
            % A setter function for featureType
            
            if ~ischar(type)
                error('Feature type must be characters. ');
            end
                       
            if strcmpi(type, 'mean')
                obj.featureType = obj.featureIsMean;
                return;
            end
            
            if strcmpi(type, 'variance')
                obj.featureType = obj.featureIsVariance;
                return;
            end
            
            if strcmpi(type, 'energy')
                obj.featureType = obj.featureIsEnergy;
                return;
            end
            
            if strcmpi(type, 'kurtosis')
                obj.featureType = obj.featureIsKurtosis;
                return;
            end
            
            if strcmpi(type, 'regression')
                obj.featureType = obj.featureIsRegression;
                return;
            end 
            
            error('Undefined feature type. ');
            
        end % function setFeatureType
    
    end % class utilities
    
end % class FeatureExtractionMean