
addpath(genpath('/network/rit/home/mh249156/EEG/toolbox/multichannelClassifier/'));
pathToData = 'D:\Research\Datasets\Metamers_Dataset\data\_SSVEP_EXP2';
targetsFrequencies = 10;
samplingFrequency = 256;
harmonicsCount = 3;
prallelProcessingFlag = true;
subject = 46;
% 
% [~, exp2Data] = loadData(pathToData, subject);
% exp2Data = squeeze(exp2Data(1, :, :, :, :));
[exp2Data, ~, signalsFlat] = loadData(pathToData, subject);
temp = transpose(signalsFlat{1});
flatSignal = zeros([1, size(temp)]);
flatSignal(1, :, :) = temp(:, :);
exp2Data = permute(exp2Data, [2, 1, 3, 4]);

exp2Data = reshape(exp2Data, [size(exp2Data, 1)*size(exp2Data, 2),...
    size(exp2Data, 3), size(exp2Data, 4)]);

exp2Data = permute(exp2Data, [1, 3, 2]);


sums = zeros(size(exp2Data));
bands1 = 5;
bands2 = bands1 + 40;

for i = 1:length(bands1)
% % Filter Data 
d = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',bands1(i),'HalfPowerFrequency2',bands2(i), ...
    'SampleRate',256);
s = size(exp2Data);
X = transpose(reshape(exp2Data, [s(1)*s(2), s(end)]));
filtered_data = filtfilt(d, X);
sums = sums + reshape(filtered_data', s);
end

exp2Data = sums(:, :, 35:end-25); 


% electrodes = 56:64;
% keptTimes = zeros(240, length(electrodes));
% endSample = 750;
% targetsFrequencies = [8:15, 8.2:15.2, 8.4:15.4, 8.6:15.6, 8.8:15.8];
% samplingFrequency = 250;
% harmonicsCount = 3;
% prallelProcessingFlag = true;
% 
% benchmarkData = loadDatabase('default path', 6); 
% benchmarkData = benchmarkData.data(electrodes, 1:endSample, 3, :);
% benchmarkData = reshape(benchmarkData, ...
%     [size(benchmarkData, 1),...
%     size(benchmarkData, 2),...
%     size(benchmarkData, 3) * size(benchmarkData, 4)]);
% currentSignal = permute(benchmarkData, [3, 1, 2]);
% exp2Data = currentSignal; 

% myMisc = FeatureExtractionMisc();
% myMisc = setupFeatureExtractionMethod(myMisc,...
%     'samplingFrequency', samplingFrequency,...
%     'featureType', 'mean',...
%     'parallelProcessing', prallelProcessingFlag);

% myARMA = FeatureExtractionARMA();
% myARMA = setupFeatureExtractionMethod(myARMA,...
%     'samplingFrequency', samplingFrequency,...
%     'orderAR', 8, ...
%     'orderMA', 0,...
%     'parallelProcessing', prallelProcessingFlag);
% 
% mySpectralCD = ChangeDetectorARMA();
% mySpectralCD = setupChangeDetector(mySpectralCD,...
%     'featureExtractor', myARMA);
% [~, allScores, windowTimes] = detectChange(mySpectralCD, exp2Data([1, 11, 21, 31], 11:15, 70:end-70));


% for i = 1:size(allScores, 1)
%     trials = [1, 11, 21, 31];
%     indexes = 1:size(allScores, 3);
%     indexes = indexes / 256;
%     currentSignal = allScores(i, :, :);
%     currentSignal = squeeze(currentSignal);
%     f = figure('visible', 'off');
%     plot(indexes, currentSignal', 'LineWidth', 4);
%     xlabel('time (s)');
%     ylabel('Log likelihood');
%     xticks(1:6);
%     xticklabels(1:6);
%     grid on;
%     set(gca, 'FontSize', 65);
%     set(gcf, 'Position',  [1, 1, 2160, 1440]);
%     myLeg = legend({'E11', 'E12', 'E13', 'E14', 'E15'});
%     myLeg.NumColumns = 2;
%     title(['Subject ', int2str(subject), ' - Trial ', int2str(trials(i))]);
%     saveas(f, ['S', int2str(subject), '_AR_metamers_Trial', int2str(trials(i)), '.png']);
% end % for i
% 


















% % subbands = [8, 88;...
% %     16, 88; 24, 88; 32,...
% %     88; 40, 88; 48, 88;...
% %     56, 88;];
% % 
% % myFBCCA = FeatureExtractionCCA();
% % myFBCCA = setupFeatureExtractionMethod(myFBCCA,...
% %     'samplingFrequency', samplingFrequency,...
% %     'harmonicsCount', harmonicsCount, ...
% %     'targetsFrequencies', targetsFrequencies,...
% %     'isFilterBank', true,...
% %     'subbands', subbands,...
% %     'parallelProcessing', prallelProcessingFlag);
% % 
% % myCCA = FeatureExtractionCCA();
% % myCCA = setupFeatureExtractionMethod(myCCA,...
% %     'samplingFrequency', samplingFrequency,...
% %     'harmonicsCount', harmonicsCount, ...
% %     'targetsFrequencies', targetsFrequencies,...
% %     'isFilterBank', false,...
% %     'parallelProcessing', prallelProcessingFlag);
% % 
% % myMSI = FeatureExtractionMSI();
% % myMSI = setupFeatureExtractionMethod(myMSI,...
% %     'samplingFrequency', samplingFrequency,...
% %     'harmonicsCount', harmonicsCount, ...
% %     'targetsFrequencies', targetsFrequencies,...
% %     'embeddingDimension', 1,...
% %     'delaySize', 1,...
% %     'parallelProcessing', prallelProcessingFlag);
% % 
% % myMEC = FeatureExtractionMEC();
% % myMEC = setupFeatureExtractionMethod(myMEC,...
% %     'samplingFrequency', samplingFrequency,...
% %     'harmonicsCount', harmonicsCount, ...
% %     'targetsFrequencies', targetsFrequencies,...
% %     'parallelProcessing', prallelProcessingFlag);
% % 
% % addpath(genpath('/network/rit/home/mh249156/EEG/toolbox/multichannelClassifier/'));
% % pathToData = 'D:\Research\Datasets\Metamers_Dataset\data\_SSVEP_EXP2';
% % targetsFrequencies = 10;
% % samplingFrequency = 256;
% % harmonicsCount = 3;
% % prallelProcessingFlag = true;
% % % 
% % [~, exp2Data] = loadData(pathToData, 46);
% % exp2Data = squeeze(exp2Data(1, :, :, :, :));
% 
% % [exp2Data, ~, signalsFlat] = loadData(pathToData, 46);
% % temp = transpose(signalsFlat{1});
% % flatSignal = zeros([1, size(temp)]);
% % flatSignal(1, :, :) = temp(:, :);
% % exp2Data = permute(exp2Data, [2, 1, 3, 4]);
% % 
% % exp2Data = reshape(exp2Data, [size(exp2Data, 1)*size(exp2Data, 2),...
% %     size(exp2Data, 3), size(exp2Data, 4)]);
% % 
% % exp2Data = permute(exp2Data, [1, 3, 2]);
% % 

% originalSize = size(exp2Data);
% signal = reshape(exp2Data, [size(exp2Data, 1)*size(exp2Data, 2), size(exp2Data, 3)]);
% signal = transpose(signal);
% exp2Data = transpose(resample(signal, 10, 1));
% exp2Data = reshape(exp2Data, [originalSize(1), originalSize(2), size(exp2Data, 2)]);

myMisc = FeatureExtractionMisc();
myMisc = setupFeatureExtractionMethod(myMisc,...
    'samplingFrequency', samplingFrequency,...
    'featureType', 'regression', ...
    'parallelProcessing', prallelProcessingFlag);

% myCD = ChangeDetectorGLR();
% myCD = setupChangeDetector(myCD, ...
%     'featureExtractor', myMisc,...
%     'windowType', 'sliding',...
%     'windowLength', 256,...
%     'windowOverlap', 0.9);
% 
myCDMD = ChangeDetectorGLRMD();
myCDMD = setupChangeDetector(myCDMD, ...
    'featureExtractor', {myMisc},...
    'windowType', 'sliding',...
    'windowLength', 256,...
    'windowOverlap', 0.9);

% % myCDMD = ChangeDetectorGLRMD();
% % myCDMD = setupChangeDetector(myCDMD, ...
% %     'featureExtractor', {myMSI},...
% %     'windowType', 'sliding',...
% %     'windowLength', 256,...
% %     'windowOverlap', 0.8);
% 
[~, allScores, windowTimes, features] = detectChange(myCDMD, exp2Data([1, 11, 21, 31, 41, 51], :, :));
% times = analyzeStartingTime(allScores, windowTimes);
% maxLikelihoods = getMaximumLikelihoods(allScores);
% 
% % for i = 1:3
% %     idx = [3:54; 57:108; 111:162];
% %     f = figure('visible', 'off');
% %     experimentData = times(idx(i, :));
% %     plot(experimentData, 'LineWidth', 4);
% %     xlabel('Trial Number');
% %     ylabel('Estimated Time (s)');
% %     title(['S46 - Trial ', int2str(i)]);
% %     grid on;
% %     ylim([0 6]);
% %     xlim([1 52]);
% %     xticks(1:7:length(experimentData));
% %     xticklabels(3:7:54);
% %     set(gca, 'FontSize', 65);
% %     set(gcf, 'Position',  [1, 1, 1440, 1440]);
% %     saveas(f, ['s46_timeVstrial_trial', int2str(i), '_MEC.fig']);
% %     saveas(f, ['s46_timeVstrial_trial', int2str(i), '_MEC.png']);
% % end
% 
% for i = 1:3
%     idx = [3:54; 57:108; 111:162];
%     f = figure('visible', 'off');
%     experimentData = maxLikelihoods(idx(i, :));
%     plot(experimentData, 'LineWidth', 4);
%     xlabel('Trial Number');
%     ylabel('Max Log Likelihood');
%     title(['S46 - Trial ', int2str(i)]);
%     grid on;
%     xlim([1 52]);
%     xticks(1:7:length(experimentData));
%     xticklabels(3:7:54);
%     set(gca, 'FontSize', 65);
%     set(gcf, 'Position',  [1, 1, 1440, 1440]);
%     saveas(f, ['s46_timeVsLikelihood_trial', int2str(i), '_MSI.fig']);
%     saveas(f, ['s46_timeVsLikelihood_trial', int2str(i), '_MSI.png']);
% end
% 
% close all
% 
% 
% function times = analyzeStartingTime(allScores, windowTimes)
% 
%     times = zeros(size(allScores, 1), 1);
%     
%     for i = 1:size(allScores, 1)
%         [~, maxLocation] = max(allScores(i, :));
%         times(i) = windowTimes(2, maxLocation);        
%     end  % for i
% 
% end % function analyzeStartingTime
% 
% 
% function maximumLikelihooods = getMaximumLikelihoods(allScores)
% 
%     maximumLikelihooods = zeros(size(allScores, 1), 1);
%     
%     for i = 1:size(allScores, 1)
%         [maxValue, maxLocation] = max(allScores(i, :));
%         maximumLikelihooods(i) = maxValue;
%     end % for i 
% 
% end


