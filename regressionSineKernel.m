function [score, residuals] =  regressionSineKernel(signal, fs, targetFrequency)
    t = 1:size(signal, 2);
    E = size(signal, 1);
    iterationsCount = 25;
    freqs = 10;
    score = zeros(length(freqs), 1);

    for i = 1:length(freqs)
        f = freqs(i);
        sin1 = sin(1*2*pi*f*t/fs);
        cos1 = cos(1*2*pi*f*t/fs);
        sin2 = sin(2*2*pi*f*t/fs);
        cos2 = cos(2*2*pi*f*t/fs);
        sin3 = sin(3*2*pi*f*t/fs);
        cos3 = cos(3*2*pi*f*t/fs);
        H = [sin1; cos1; sin2; cos2; sin3; cos3];
        W = estimateW(signal, H, iterationsCount);
        X = estimateX(signal, H, W);
        residuals = W - X;
        score(i) = computeScore(signal, W, X, H);    
    end
end


function W = estimateW(s, H, iterationsCount)
    E = size(s, 1);
    W = ones(E, 1)/E;
    for iteration = 1:iterationsCount

        for j = 1:E
            k1 = zeros(size(H, 1));
            k2 = zeros(size(H, 1), 1);

            for i = 1:size(s, 2)
                k1 = k1 + H(:, i) * H(:, i)';
                k2 = k2 + s(j, i) * H(:, i);        
            end % for i

            k1 = inv(k1);
            k0 = k1 * k2 / (E * sum(s(j, :)));

            k3 = 0;
            for m = 1:E
                for i = 1:size(s, 2)
                    k3 = W(m) * s(m, i) * H(:, i)';
                end % for i
            end % for j

            newW(j) = k3 * k0;

        end % for j

        W = newW / sum(newW);
    end % for iteration

end % function estimateW

function X = estimateX(s, H, W)
    E = size(s, 1);
    X = zeros(size(H, 1));
    
    for i = 1:size(s, 2)
        X = X + H(:, i) .* H(:, i)';    
    end % for i
    
    X = inv(X)/E;
    k = zeros(1, size(H, 1));

    for i = 1:size(s, 2)
        for j = 1:E
            k = k + W(j) * s(j, i) * H(:, i)';
        end
    end

    X = k * X;
    X = X';
end % function estimateX

function score = computeScore(s, W, X, H)
    score = 0;
    for j = 1:size(s, 1)
        for i = 1:size(s, 2)
            score = score + (W(j) * s(j, i) * X'*H(:, i))^2;
        end 
    end 
end % function computeScore