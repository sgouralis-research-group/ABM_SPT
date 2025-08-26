function [C1, C2] = cfun(t, tstar, T, D)
    % function now returns two blocks C1, C2. The entire covariance matrix is
    % C1 0
    % 0 C2
    % in the rare case that we want just the Lambda matrix without
    % multiplying by 2D, we can call the function without a D argument
    if ~exist('D', 'var')
        C1 = min(abs(t(t<T)-T), abs(tstar(tstar<T)-T));
        C2 = min(abs(t(t>=T)-T), abs(tstar(tstar>=T)-T));
    else
        C1 = 2*D*min(abs(t(t<T)-T), abs(tstar(tstar<T)-T));
        C2 = 2*D*min(abs(t(t>=T)-T), abs(tstar(tstar>=T)-T));
    end
