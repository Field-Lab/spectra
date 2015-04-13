function [sigma, mu] = getNoise(data)
% Matricialized version of noise calculation
% data should be 1 column per electrode, size(data) = nSamples,nElectrodes

data = double(data); % Cast data - for compatibility calls from uint8
% TODO Consider going to single for speed
data2 = data.^2; % Pre computing squared data only once

% Extracting parameters
nSamples = size(data,1);
nElectrodes= size(data,2);

% Output arrays
sigma = inf(1,nElectrodes);
mu = zeros(1,nElectrodes);

oneMat = ones(nSamples,1);

for el = 1:nElectrodes % Working individually on electrodes - the number of while iterations required can be very different
    el; % Speed debug
    % Electrode setup - Reset counters
    hitCount = 0;
    sumX = 0;
    sumX2 = 0;
    
    while true % Iterating for noise rms value convergence
        
        mask = abs(data(:,el) - mu(el)*oneMat) < 3*sigma(el)*oneMat; % Masking outliers out
        % Update counters
        hitCount = hitCount + sum(mask);
        sumX = sumX + sum(data(:,el).*mask);
        sumX2 = sumX2 + sum(data2(:,el).*mask);
        
        % Update average and std
        sigmaNew = sqrt((sumX2 - sumX.^2./hitCount)./(hitCount-1));
        mu(el) = sumX./hitCount;
        
        % Check for convergence
        if abs(sigma(el)- sigmaNew) < 0.01 % WARNING HARDOCED 0.01
            sigma(el) = sigmaNew;
            break
        else
            sigma(el) = sigmaNew;
        end
    end % while
end % el

end

