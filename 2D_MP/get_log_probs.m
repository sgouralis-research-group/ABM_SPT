function log_P = get_log_probs(v,D,T_vec,X_vec,Y_vec,x_cell,y_cell,params)
% log_P = log( post,like,priors )
if nargin<1
    log_P = 9;
    return
end

log_P = zeros(get_log_probs,1);

for i = 1:params.M
    %% Likelihood
    % solve the 2d system by hand using the SVD of H
    x = params.wx_cell{i}-x_cell{i};
    y = (params.wy_cell{i}-y_cell{i} - params.h_21_U_cell{i}.*x_cell{i}./params.h_11_U_cell{i}) ...
        ./ (params.h_22_U_cell{i} - params.h_21_U_cell{i}.*params.h_12_U_cell{i}./params.h_11_U_cell{i});
    
    x = (x - params.h_12_U_cell{i}.*y)./params.h_11_U_cell{i};

    log_P(2) = sum(  x.^2 ./ params.h_1_S_cell{i} ...
                   + y.^2 ./ params.h_2_S_cell{i}) + log_P(2);
  
    % think this could be vectorized outside the loop
    log_P(2) = log_P(2) - params.N_vec(i)*log(v); 

    [c1, c2] = cfun(params.t_cell{i}, params.t_cell{i}', T_vec(i), D);

    eigs = eig(blkdiag(c1,c2));
    
    L = chol(blkdiag(c1,c2), 'lower');

    L_inv_xX = forwardSubstitution(L, x_cell{i}-X_vec(i));
    L_inv_yY = forwardSubstitution(L, y_cell{i}-Y_vec(i));

    % x
    log_P(7) =  sum(L_inv_xX.^2)  ...
                - 0.5 * sum(log(eigs)) ...
                + log_P(7);

    % y
    log_P(9) =  sum(L_inv_yY.^2) ...
                -0.5 * sum(log(eigs)) ...
               +log_P(9);
end
%v
log_P(3) = log((((params.v_prior_A-1)*params.v_prior_ref)^params.v_prior_A)...
                /gamma(params.v_prior_A)*v^(-(params.v_prior_A-1))...
                *exp(-(params.v_prior_A-1)*params.v_prior_ref/v));

% D
log_P(4) = log((((params.D_prior_A-1)*params.D_prior_ref)^params.D_prior_A)...
                /gamma(params.D_prior_A)*D^(-(params.D_prior_A-1))...
                *exp(-(params.D_prior_A-1)*params.D_prior_ref/D));

% T
log_P(5) = - sum ( log(params.T_prior_max_vec-params.T_prior_min_vec)) ;

% X
log_P(6) = -1/2 * sum ( log(params.X_prior_U_vec).*((-1/2)*(X_vec-params.X_prior_ref_vec).^2./params.X_prior_U_vec) );

% Y
log_P(8) = -1/2 * sum ( log(params.Y_prior_U_vec).*((-1/2)*(Y_vec-params.Y_prior_ref_vec).^2./params.Y_prior_U_vec) );

%% Posterior
log_P(1) = sum(log_P(2:end));




