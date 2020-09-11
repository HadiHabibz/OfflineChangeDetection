classdef ChangeDetectorARMA < ChangeDetector
    
  properties(Access = protected)
        % Empty
        reference;
    end % class protected data members
    
    methods(Access = public)
        
        function obj = ChangeDetectorARMA()
            % Constructor
            
            obj@ChangeDetector();   
        end % class constructor
        
        function [obj, allScores, windowsTime] = detectChange(obj, signal)
           
            windowsTime = 0;
            tic
            allScores = zeros(size(signal));
            
            for i = 1:size(signal, 1)
                
                obj.featureSignal = signal(i, :, :);

%                 d = designfilt('bandpassiir','FilterOrder',10, ...
%                     'HalfPowerFrequency1',8,'HalfPowerFrequency2',32, ...
%                     'SampleRate',256);
%                 obj.featureSignal = filter(d, obj.featureSignal);

                % Test all hypotheses and compute the log likelihood for each
                % of them.
                [allScores(i, :, :),...
                    windows] = computeAllPossibleChangeTimes(obj);  
                
            end % for i
            toc
                       
        end % function detectChange

    end % class public services
    
    methods(Access = protected)

        function [allScores, windows] = ...
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
            allScores = zeros(size(obj.featureSignal, 2), k);
            
            if k <= 1
                error('process too short.\n');
            end
            
            maxOrder = max(getOrderAR(obj.featureExtractor),...
                getOrderMA(obj.featureExtractor));
            
            
            electrodesCount = size(obj.featureSignal, 2);
            for i = 1:electrodesCount
                Y = obj.featureSignal(1, i, :);
                Y = squeeze(Y);                
                L = 10;                
                p = getOrderAR(obj.featureExtractor);              
                [eS, ES, rS, RpS, CpS, dS] = ladderAlgorithm(obj, Y', L, p);
                
                for t = L+1:length(Y)
                    L_t = t;
                    [eG, EG] = ladderAlgorithmGrowing(obj,...
                        Y', L, p, L_t, rS, ES, RpS, CpS, dS);
                    Hst = L * log( ES(p, t) / L );
                    Hgt = L_t * log( EG(p, t) / L_t );
                    Hgtl = (L_t - L) * log( EG(p, t-L) / (L_t - L) );
                    d(t) = Hgt - Hst - Hgtl;
                    
                end % for t
            end % for i
            windows = 0;
            
%             for j = 125+maxOrder:k-maxOrder-125
%                    
%                 score = computeScore(obj, j, k);
%                 allScores(:, j) = score;
%                 windows(j) = j;
%                 
%             end % for j

            
        end % end function computeAllPossibleChangeTimes
        
        function distance = computeScore(obj, j, k)
            
            distance = zeros(size(obj.featureSignal, 2), 1);
            AROrder = getOrderAR(obj.featureExtractor);
            s = j-125;

            
            for i = 1:size(obj.featureSignal, 2)   
                 
                referenceWindow = ...
                    obj.featureSignal(:, i, s:j-1);
                referenceWindow = squeeze(referenceWindow);
                referenceWindow = referenceWindow - mean(referenceWindow);
                [arReference, err] = aryule(referenceWindow, AROrder);               
                
                r = j + 125;
                       
                testWindow = ...
                    squeeze(obj.featureSignal(:, i, j:r));
                testWindow = testWindow - mean(testWindow);
                [arTest, ert] = aryule(testWindow, AROrder);
              
                
                pooledWindow = ...
                    obj.featureSignal(:, i, s:r);
                pooledWindow = squeeze(pooledWindow);
                pooledWindow = pooledWindow - mean(pooledWindow);
                [arPooled, erp] = aryule(pooledWindow, AROrder);                 
                
%                 testErrorR = filter(arReference, 1, referenceWindow) + referenceWindow;
%                 testErrorT = filter(arTest, 1, testWindow) + testWindow; 
%                 pooledError = filter(arPooled, 1, pooledWindow) + pooledWindow;
                
%                 erp = sum(pooledError.^2);
%                 err = sum(testErrorR.^2);
%                 ert = sum(testErrorT.^2);
                distance(i) = length(pooledWindow) * log(erp/length(pooledWindow)) -...
                    length(referenceWindow) * log(err/length(referenceWindow)) - ...
                    length(testWindow) * log(ert/length(testWindow));      

%                 distance(i) = err / ert;
                
            end % for i
         
        end % function computeScore            
        
        function [e, E,  r, Rp, Cp, d] = ladderAlgorithm(obj, Y, L, p)
            % Ensure Y is a row vector
            tmax = length(Y)+1;
            r = ones(p, tmax);
            q = ones(p, tmax);
            Ep = ones(p, tmax);
            Rp = ones(p, tmax);
            Cp = ones(p, tmax);
            sigma = ones(p, tmax);
            pi = ones(p, tmax);
            E = ones(p, tmax);
            e = ones(p, tmax);
            d = ones(p, tmax);
            R = ones(p, tmax);
            C = ones(p, tmax);
            y = [zeros(1, L), Y];
            sigma(:, 1) = 1;
            pi(:, 1) = 1;
            
            for t = 1:length(Y)
                past_t = t;
                current_t = t+1;                
                e(1, current_t) = Y(t);
                r(1, current_t) = Y(t);
                d(1, current_t) = y(t);
                q(1, current_t) = y(t);
                sigma(1, current_t) = 1;
                pi(1, current_t) = 1;
                for n = 1:p
                    temp = Ep(n, past_t) + (e(n, current_t)^2)/sigma(n, past_t);
                    E(n, current_t) = temp;
                    
                    temp = E(n, current_t) - (d(n, current_t)^2)/pi(n, past_t);
                    Ep(n, current_t) = temp;
                    
                    temp = Rp(n, past_t) + (r(n, past_t)^2)/sigma(n, past_t);
                    R(n, current_t) = temp;
                    
                    temp = R(n, current_t) - (q(n, past_t)^2)/pi(n, past_t);
                    Rp(n, current_t) = temp;
                    
                    temp = Cp(n, past_t) + (e(n, current_t)*r(n, past_t))/sigma(n, past_t);
                    C(n, current_t) = temp; 
                    
                    temp = C(n, current_t) - (d(n, current_t)*q(n, past_t))/pi(n, past_t);
                    Cp(n, current_t) = temp;
                                        
                    temp = e(n, current_t) - (r(n, past_t)*C(n, current_t))/R(n, current_t);
                    e(n+1, current_t) = temp;
                    
                    temp = d(n, current_t) - (q(n, past_t)*C(n, current_t))/R(n, current_t);
                    d(n+1, current_t) = temp;                    
                    
                    temp = r(n, past_t) - (e(n, current_t)*C(n, current_t))/E(n, current_t);
                    r(n+1, current_t) = temp;
                    
                    temp = q(n, past_t) - (d(n, current_t)*C(n, current_t))/E(n, current_t);
                    q(n+1, current_t) = temp;
                                        
                    temp = sigma(n, past_t) - (e(n, current_t)^2)/E(n, current_t);
                    sigma(n+1, current_t) = temp;
                    
                    temp = pi(n, past_t) - (d(n, current_t)^2)/(E(n, current_t));                    
                    pi(n+1, current_t) = temp;
                    
                end % for n
            end % for t
            
%             E = E(:, L+1:end);
%             e = e(:, L+1:end);
            
        end % function ladderAlgorithm
%         
%         function [e, E] = ladderAlgorithmGrowing(obj, y, L, p, L_t, r_sliding, E_sliding, Rp_sliding, Cp_sliding, d_sliding)
%             tmax = length(y);
%             r = zeros(p, tmax);
%             e = zeros(p, tmax);
%             sigma = zeros(p, tmax);
%             pi = zeros(p, tmax);
%             h = zeros(p, tmax);
%             E = zeros(p, tmax);
%             d = zeros(p, tmax);
%             ep = zeros(p, tmax);
%             Ep = zeros(p, tmax);
%             sigmap = zeros(p, tmax);
%             Cp = zeros(p, tmax);
%             Rp = zeros(p, tmax);
%             
%             for t = 2:tmax
%                 
%                 if L_t <= L
%                     return;
%                 end 
%                 
%                 if L_t <= L + p
%                     n = L_t - L;
%                     
%                     if n == 1                       
%                         r(1, t-1) = r_sliding(1, t-1);
%                         sigma(1, t-1) = 1;
%                         pi(1, t-1) = 1;
%                         h(1, t-1) = 0;                                                                       
%                     end % if n
%                     
%                     E(n, t-1) = E_sliding(n, t-1);
%                     Rp(n, t-1) = Rp_sliding(n, t-1);
%                     Cp(n, t-1) = Cp_sliding(n, t-1);
%                     d(n, t-1) = d_sliding(n, t-1);
%                 end
%                 
%                 e(1, t) = y(t);
%                 r(1, t) = y(t);
%                 sigma(1, t) = 1;
%                 pi(1, t) = 1;
%                 h(1, t) = 0;
%                 
%                 for n = 1:min(p, L_t - L)
%                     E(n, t) = E(n, t-1) + e(n, t)^2 / sigma(n, t-1);
%                     d(n, t) = d(n, t-1) + (e(n, t)*h(n, t-1)) / sigma(n, t-1);
%                     ep(n, t) = e(n, t) - (d(n, t)*h(n, t-1)) / pi(n, t-1);                    
%                     Ep(n, t) = E(n, t) - d(n, t)^2 / pi(n, t-1);
%                     sigmap(n, t) = sigma(n, t-1) - h(n, t-1)^2 / pi(n, t-1);
%                     Rp(n, t) = Rp(n, t-1) + r(n, t-1)^2 / sigmap(n, t);
%                     Cp(n, t) = Cp(n, t-1) + ep(n, t)*r(n, t-1)/sigmap(n, t);
%                     e(n+1, t) = ep(n, t) - r(n, t-1)*Cp(n, t)/Rp(n, t);
%                     r(n+1, t) = r(n, t-1) - ep(n, t)*Cp(n, t)/Ep(n, t);
%                     sigma(n+1, t) = sigma(n, t-1) - e(n, t)^2 / E(n, t);
%                     pi(n+1, t) = pi(n, t-1) - d(n, t)^2 / E(n, t);
%                     h(n+1, t) = h(n, t-1) - d(n, t)*e(n, t) / E(n, t);                    
%                 end % for n
%                 
%       
%                 
%             end % for t
%             
%         end % function ladderAlgorithmGrowing
        
    end % class utilities   
    
end % class ChangeDetectorARMA