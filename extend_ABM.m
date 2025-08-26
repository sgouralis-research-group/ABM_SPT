function [x_ext,t_ext] = extend_ABM(t_add, x, t, D, X, T)
    [t_ext,I] = sort([t', t_add']);
    t_ext = t_ext';
   
    [C1tt, C2tt] = cfun(t', t, T, D);             % sig22
    [C1tta, C2tta] = cfun(t, t_add', T, D);       % sig21
    [C1tat, C2tat] = cfun(t_add, t', T, D);       % sig12
    [C1tata, C2tata] = cfun(t_add', t_add, T, D); % sig11

    opts.SYM=true; opts.POSDEF=true;
        
    % (x1 | x2 = a) ~ N(MuBar, SigmaBar)
    % MuBar = mu_1 + sig12*sig22^-1*(a-mu_1)
    % SigBar = Sig11 - Sig12*Sig22^-1*Sig21 
    MuBar = X + blkdiag(C1tat, C2tat)*linsolve(blkdiag(C1tt, C2tt), x-X, opts);
    SigBar = blkdiag(C1tata, C2tata) - blkdiag(C1tat, C2tat)* ...
             linsolve(blkdiag(C1tt, C2tt), blkdiag(C1tta, C2tta), opts);

    [U, S, ~] = svd(SigBar);
    L = U*sqrt(S);
    
    % combine old x with the new sample x and sort them correctly
    x_ext = [x', (MuBar+L*randn(length(t_add), 1))']';
    x_ext = x_ext(I);