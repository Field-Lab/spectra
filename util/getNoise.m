function [sigma, mu] = getNoise(data)
    %GETNOISE Computes noise floor by a discard-and-recurse process
    %
    %   Input:
    %       data: nElectrodes x nSamples double array
    %
    % Author -- Vincent Deo -- Stanford University -- August 27, 2015
    
    % NECESSARY CAST - otherwise in single - precision overflow in covariance summing - neg sqrt...
    if ~isa(data,'double')
        data = double(data);
    end
    
    validateattributes(data,{'double'},{'2d','nonempty','nonsparse','real','finite','nonnan'},'','data',1);
    if size(data,2) <= 1
        throw(MException('','getNoise:cannot compute noise floor on only one sample'));
    end
    
    data2 = data.^2; % Pre computing squared data only once
    
    % Extracting parameters
    nElectrodes= size(data,1);
    
    % Output arrays
    sigma = inf(nElectrodes,1);
    mu = zeros(nElectrodes,1);
    
    for el = 1:nElectrodes % Working individually on electrodes - the number of while iterations required can be very different
        
        % Electrode setup - Reset counters
        hitCount = 0;
        sumX = 0;
        sumX2 = 0;
        
        dataEl = data(el,:);
        data2El = data2(el,:);
        
        while true % Iterating for noise rms value convergence
            
            mask = abs(dataEl - mu(el)) < 3*sigma(el); % finding outliers
            % Chop off outliers
            dataEl(~mask) = [];
            data2El(~mask) = [];
            
            % Update counters
            hitCount = hitCount + numel(dataEl);
            sumX = sumX + sum(dataEl);
            sumX2 = sumX2 + sum(data2El);
            
            % Update average and std
            sigmaNew = sqrt((sumX2 - sumX.^2./hitCount)./(hitCount-1));
            mu(el) = sumX./hitCount;
            
            % Check for convergence
            if abs(sigma(el)- sigmaNew) < 0.01 % Hardcoded convergence compliance
                sigma(el) = sigmaNew;
                break
            else
                sigma(el) = sigmaNew;
            end
        end % while
    end % el
end