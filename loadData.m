function [data, dataSynthetic, signalsFalt] = loadData(pathToData, subjectID)
% For each participant, there are at least three experiments with 54
% trials. Each trial is six seconds with an additional 1 second
% interstimulation interval. Thus, the output data is a three dimensional
% array with size [3, 54, N, 16], where N is the number of samples and I
% assume there are 16 electrudes. Each extracted signal is about seven
% secconds. The first second is pre-stimulation and the next six seconds
% are post-stimulation. The second output, dataSynthetic, is a synthetic
% data. Each signal is 12 seconds. The first 6 seconds are with
% no stimulation and the second 6 seconds are with stimulation.
% 'dataSynthetic' has dimensions (2, 52, 2*M, 16), where M is the number of
% samples and 16 is the number of channels. When the first (second) 
% entry of the array is one (i.e., dataSynthetic(1,  m,  :, :) (two), 
% the signal is constructed as follows: The first six seconds of the signal 
% are the no-stimulation signal of the first (second) trial. The next six 
% seconds are the stimulation period of the (m+2)th tria;. 

% Construct the path up to the subject directory
if subjectID <= 9
    subjectStringID = ['S0', int2str(subjectID)];
else
    subjectStringID = ['S', int2str(subjectID)];
end % if subjectID

pathToData = fullfile(pathToData, subjectStringID);

% Obtain how many dat files are in the folder
folder = dir(fullfile(pathToData, '*.dat'));
filesCount = size(folder, 1);
experimentCounter = 1;
signalsFalt = {};

for i = 1:filesCount
    
    filename = [subjectStringID, 'S001R0', int2str(i), '.dat'];
    
    [signal,states,parameters] =...
        load_bcidat(fullfile(pathToData,filename));
    
    if parameters.SamplingRate.NumericValue ~= 256
        fprintf("Warning! Odd sampling rate for %s\n", filename);
    end
    
    [trialsCount, edges] = extractTrialsCount(states);
    
    if trialsCount ~= 54
        continue;
    end 
    
    [signals, syntheticSignals, boxes] = tokenize(signal, edges, trialsCount);
    signalsFalt(experimentCounter) = {signal};
    
    if experimentCounter == 1
        data = zeros([3, size(signals)]);
        dataSynthetic = zeros([3, size(syntheticSignals)]);
    end
    
    data(experimentCounter, :) = signals(:);
    dataSynthetic(experimentCounter, :) = syntheticSignals(:);
    experimentCounter = experimentCounter + 1;
    
    if experimentCounter > 3
        break;
    end
  
% Use these to see if the extracted signals match captured
% data well enough!
%     figure;
%     index = 1:length(signal);
%     plot(index, states.DigitalInput1, 'b', index, boxes, 'g');
%     ylim([0, 2]);
    
end % for i
    
end % function loadData

function [tokens, syntheticTokens, boxes] = tokenize(signal, edges, trialsCount)
    
   samplingFreq = 256;
   
   % The stimulation period is not exactly six seconds. It is a few milli
   % seconds off. Likewise, the pre-stimulation interval is not exactly one
   % second. This variable disregards this many samples from the signal to
   % make sure there is no overlap between extracted epochs. 
   grace = 36;
   
   samplesCount = 7 * samplingFreq - grace;
   samplesCountNoPreStim = 6 * samplingFreq - grace/2;
   tokens = zeros(trialsCount, samplesCount, size(signal, 2));
   syntheticTokens = ...
       zeros(2, trialsCount-2, 2*samplesCountNoPreStim, size(signal, 2));
   startSamples = find(edges == 1);
   lastIndexPrevious = 0;
   boxes = zeros(size(signal));
   
   for i = 1:trialsCount
       % Include one second (ish) of pre-stimulation
       firstIndex = startSamples(i) - (samplingFreq-(grace/2));
       firstIndexNoPreStim = startSamples(i);
       
       % Include 6 seconds of stimulation
       lastIndex = startSamples(i) + (6*samplingFreq - (grace/2) - 1);
       
       if firstIndex <= lastIndexPrevious
           fprintf('Warning! Overlapping tokens.\n');
       end
       
       if firstIndex <= 0 || lastIndex >= length(signal)
           fprintf('Warning! Tokens exceed signal boundaries!\n');
           continue;
       end 
       
       lastIndexPrevious = lastIndex;
       boxes(firstIndex:lastIndex) = 1;       
       tokens(i, :, :) = signal(firstIndex:lastIndex, :);
       
       % The first two trials are baselines
       if i == 1
           noStimulationSourceOff =...
               signal(firstIndexNoPreStim:lastIndex, :);
       elseif i == 2
           noStimulationSourceOn = ...
               signal(firstIndexNoPreStim:lastIndex, :);
       else
           syntheticTokens(1, i-2, :, :) = ...
               [noStimulationSourceOff;...
               signal(firstIndexNoPreStim:lastIndex, :)];
           syntheticTokens(2, i-2, :, :) = ...
               [noStimulationSourceOn;...
               signal(firstIndexNoPreStim:lastIndex, :)];
       end 
             
   end % for i
  
end % function tokenize

function [trialsCount, edges] = extractTrialsCount(states)
    signal = transpose(states.DigitalInput1);
    signalShifted = [signal(end), signal(1:end-1)];
    edges = signal - signalShifted;
    edges = (edges == 1);
    trialsCount = sum(edges);    
end % function extractTrialsCount

% /network/rit/home/mh249156/EEG/metamersDataset/dataset/metamers_data/data/EXP1_SSVEP