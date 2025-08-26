function x = sample_ABM(t, D, X, T)
    [C1, C2] = cfun(t, t', T, D);

    % svd of each block
    [U1, S1, ~] = svd(C1, 'vector');
    [U2, S2, ~] = svd(C2, 'vector');
    
    % use implicit expansion to multiply row vectors S1', S2' by matrix U1, U2 resp.
    x = [(X+(sqrt(S1)'.*U1)*randn(length(U1), 1))' , (X+(sqrt(S2)'.*U2)*randn(length(U2), 1))']';