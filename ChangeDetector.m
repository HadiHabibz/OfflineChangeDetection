classdef ChangeDetector
    % A generic class for all change detection algorithms
                
    properties(Access = protected)
        
        % Must be of type MultichannelClassifier
        featureExtractor;
        
        % The number of signals. If there is only one signal (e.g., in
        % online experiments), this is set to 1. Otherwise, for batch
        % processing, this must be set to a value other than one
        signalsCount;
        
        % The number of samples in each window.
        % The percentage of overlap between windows.
        % Increasing the number of samples in a window helps by reducing
        % the noice. For example, we use CCA, having more samples increases
        % the performance of CCA. However, if we use big windows that
        % include many samples, the number of features that we can extract
        % reduces. Thus there is an overlap. Likewise, increasing
        % windowOverlap helps because it generates mroe features but the
        % generated features will be more dependent. In our analysis, we
        % assume that samples are independent. 
        windowLength;
        windowOverlap;
        
        % The number of samples in time domain
        samplesCount;
        
        % The number of features in the feature signal
        featuresCount;
        
        % The signal containing the features of all signals.
        % This will have size NS by W, where NS is the number of signals to
        % be processed times the number of targets and W is the length of
        % the extracted feature signal that depends on the length of the
        % signal in time domain, window size, and overlap. 
        featureSignal;
               
        % The sampling rate of data
        samplingFrequency;
        
        % The dimension of the features matrix
        featuresMatrixDim;
        
        % Determines the type of the window.
        % 'sliding' uses a moving window with a given size and overlap to
        % tokenize the signal.
        % 'extending' extends the window length by a given number of
        % samples. The start of the window remains the same. 
        windowType;
        slidingWindow;
        extendingWindow;
        
    end % class private data members
    
    methods(Access = public)
         
        function obj = ChangeDetector()
            % Constructor
            % empty
            obj.slidingWindow = 100;
            obj.extendingWindow = 101; 
        end % Class Constructor
        
        function obj = setupChangeDetector(obj, varargin)
            % Parse through input parameters
            
            obj = parseMandatoryVariables(obj, varargin);
            
        end % function setupChangeDetector
        
    end % class public services
    
    methods(Access = protected)
               
        function [obj, windowsTime] = ...
                extractFeatureSignal(obj, signal)
            % Receive a time series, truncate using multiple possibly
            % overlapping windows, and apply the designated feature
            % extraction method to each window to extract features. Thus,
            % effectively, convert the time domain signal to feature space.                              
            
            % Give the first and last index of each window. 
            [startIndex, stopIndex] = getWindowIndexes(obj);            
            obj.featuresCount = length(stopIndex);      
            
            % This helps us with plotting things
            % Keep track of when in time the windows starts
            % It has three dimensions for the begining, middle, and end
            % point of the window in time. 
            windowsTime = zeros(3, obj.featuresCount);
            
            for i = 1:obj.featuresCount
                currentSignal = signal(:, :, startIndex(i):stopIndex(i));
                windowsTime(:, i) = [startIndex(i),...
                    (startIndex(i) + stopIndex(i))/2, stopIndex(i)];
                testvalues = extractFeatures(...
                    obj.featureExtractor, currentSignal);
                
                if i == 1
                    obj.featuresMatrixDim = size(testvalues);
                    obj.featureSignal = zeros(...
                        obj.signalsCount * obj.featuresMatrixDim(2),...
                        obj.featuresCount);
                end
                
                % Save extracted features of this current window 
                obj.featureSignal(:, i) = testvalues(:);
                
            end % for i
                           
            windowsTime = windowsTime ./ obj.samplingFrequency;       
            
        end % function extractFeatureSignal
        
        function [firsts, lasts] = getWindowIndexes(obj)
            
            firstIndex = 1;
            lastIndex = obj.windowLength;
            firsts = [];
            lasts = [];
            
            overlap = obj.windowOverlap * obj.windowLength;
            overlap = floor(overlap);
                    
            while true
                
                if obj.windowType == obj.slidingWindow
                    lastIndex = firstIndex + obj.windowLength - 1; 
                end
                
                if lastIndex > obj.samplesCount
                    break;
                end
                
                firsts = [firsts, firstIndex];
                lasts = [lasts, lastIndex];
                                
                if obj.windowType == obj.slidingWindow                    
                    firstIndex = lastIndex + 1 - overlap;
                    
                elseif obj.windowType == obj.extendingWindow
                    lastIndex = lastIndex + overlap;                    
                end 
                
            end % while
            
        end % function getWindowIndexes
        
                       
        function [obj, timeSeries] = ...
                checkTimeSeriesInput(obj, newTimeSeries)               
            % Make sure the input signal has the right dimensions
            % Ideally, it must be a 3D array with dimensins (N, E, T),
            % where N is the number of signal that are to be processed in
            % batched (1 if there is a single signal), E is the number of
            % channels, and T is the number of samples.
            
            % Extract the number of signals if the dimensions are correct
            if ndims(newTimeSeries) == 3
                
                if size(newTimeSeries, 3) ~= 1
                    obj.signalsCount = size(newTimeSeries, 1);
                    timeSeries = newTimeSeries;
                    return;
                end
                
                error('Error: Samples count must exceed 1.\n');                
            end
                
            % if 2D, add a dummy dimension to make it 3D. 
            if ismatrix(newTimeSeries)
                
                if size(newTimeSeries, 2) ~= 1
                    timeSeries(1, :, :) = newTimeSeries;
                    obj.signalsCount = 1;
                    return;
                end
                
                error('Error: Samples count must exceed 1.\n');                
            end 
            
            error('Error: Times series not formatted properly. ');           
        end % end function checkTimeSeriesInput        
        
        function obj = parseMandatoryVariables(obj, options)
            % Parse the input and set variables accordingly.
            % This only parses mandatory variables
            % The following parameters must be provided:
            % 'featureExtractor' is the type of feature extraction method 
            % used to obtain the features. It must an object of class
            % MultichannelClassifier. The feature extractor must be setup
            % properly before passing to the class.
            % 'windowLength' is the number of samples in each window.
            % 'windowOverlap' is the percentage of overlap between two 
            % windows.
            % 'windowType' determines what kind of window to be used. Valid
            % inputs are 'extending' or 'sliding';
            parser = inputParser;
            parser.KeepUnmatched = true;
            addParameter(parser, 'featureExtractor', 0);
            addParameter(parser, 'windowLength', 2);
            addParameter(parser, 'windowOverlap', 0);
            addParameter(parser, 'windowType', 'sliding');
            parse(parser, options{:});
            obj.featureExtractor = parser.Results.featureExtractor;                        
            obj = setWindowLength(obj, parser.Results.windowLength);            
            obj = setWindowOverlap(obj, parser.Results.windowOverlap);     
            obj = setWindowType(obj, parser.Results.windowType);
        end % function parseMandatoryVariables
        
        function obj = setFeatureExtractor(obj, featureExtractor)
            obj.featureExtractor = featureExtractor;
        end % setter function for featureExtractor
        
        function obj = setWindowType(obj, type)
            
            if ~ischar(type)
                error('Window type must be a string.\n');
            end
            
            if strcmpi('extending', type)
                obj.windowType = obj.extendingWindow;
                return;
            end
            
            if strcmpi('sliding', type)
                obj.windowType = obj.slidingWindow;
                return;
            end
            
            error('Invalid window type. ');            
            
        end % function setWindowType
        
        function obj = setWindowLength(obj, length)            
            defaultVal = 10;
            
            if length <= 1
                fprintf('Window length too short'); 
                fprintf('Set to %d\n instead.\n', defaultVal);
                obj.windowLength = defaultVal;
                return;
            end  
            
            obj.windowLength = length;
        end % setter function for windowLength
        
        function obj = setWindowOverlap(obj, overlap)
            
            if overlap >= 0 && overlap <= 1
                obj.windowOverlap = overlap;  
                
            elseif overlap > 1 && overlap <= 100
                obj.windowOverlap = overlap / 100;                
                
            else
                fprintf('Invalid window overlap. ');
                fprintf('Set to zero by default instead.\n');
                obj.windowOverlap = 0;
            end
            
        end % setter function for windowOverlap
       
        function [obj, signal] = setOtherParameters(obj, signal)
            % After receiving the feature extractor method
            % set the remaining parameters
            
            % Ensure the input signal is 3D.
            % If 2D make it 3D by adding a leading singleton dimension.
            [obj, signal] = checkTimeSeriesInput(obj, signal);
                       
            if iscell(obj.featureExtractor)
                aFeatureExtractor = obj.featureExtractor{1};
                
            elseif length(obj.featureExtractor) == 1
                aFeatureExtractor = obj.featureExtractor;
                
            else
                error('Invalid feature extractor method.');     
                
            end % if...elseif...else
            
            obj.samplesCount = size(signal, 3);
            
            obj.samplingFrequency = ...
                getSamplingFrequency(aFeatureExtractor);
            
        end % function setOtherParameters
        
        function obj = setSamplesCount(obj, count)
            obj.samplesCount = count;
        end % setter function for samplesCount
        
    end % class protected utilities
    
end % class ChangeDetectorGLR