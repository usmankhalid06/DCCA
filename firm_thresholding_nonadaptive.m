function x_firm = firm_thresholding_nonadaptive(x, lambda1, lambda2)
    % Firm thresholding function
    % Inputs:
    %   x       : Input data (can be a vector or matrix)
    %   lambda1 : Lower threshold (soft thresholding threshold)
    %   lambda2 : Upper threshold (hard thresholding threshold)
    % Output:
    %   x_firm  : Output data after applying firm thresholding

    % Initialize the output array
    x_firm = zeros(size(x));
    
    % Apply firm thresholding
    for i = 1:length(x)
        if abs(x(i)) <= lambda1
            % Set to zero if below lower threshold (like hard thresholding)
            x_firm(i) = 0;
        elseif abs(x(i)) > lambda1 && abs(x(i)) <= lambda2
            % Apply partial shrinkage between the two thresholds
            x_firm(i) = sign(x(i)) * (lambda2 * (abs(x(i)) - lambda1) / (lambda2 - lambda1));
        else
            % Keep the value if above the upper threshold
            x_firm(i) = x(i);
        end
    end
end
