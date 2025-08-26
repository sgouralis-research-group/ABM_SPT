function sample = chainer_init_sample(params, opts)

%% counter
sample.i = 0;

%% model variables
sample.v = (params.v_prior_A-1)*params.v_prior_ref / randg(params.v_prior_A);
%sample.v = opts.ground.v;
sample.D = (params.D_prior_A-1)*params.D_prior_ref / randg(params.D_prior_A);
%sample.D = opts.ground.D;
  
sample.T_vec = params.T_prior_min_vec ...
    + (params.T_prior_max_vec-params.T_prior_min_vec) ...
    .* rand(opts.M,1);

%make sure none of the Ts are within the tolerance of timepoints
for m=1:opts.M
    while min(abs(sample.T_vec(m)-params.t_cell{m})) < params.h_tol
     sample.T_vec(m) = params.T_prior_min_vec(m) ...
        + (params.T_prior_max_vec(m)-params.T_prior_min_vec(m)) ...
        * rand;
    end
end
% sample.T_vec = opts.ground.T_vec;

sample.X_vec = params.X_prior_ref_vec ...
    + sqrt(params.X_prior_U_vec) ...
    .* randn(opts.M,1);
% sample.X_vec = opts.ground.X_vec;

sample.Y_vec = params.Y_prior_ref_vec ...
    + sqrt(params.Y_prior_U_vec) ...
    .* randn(opts.M,1);
% sample.Y_vec = opts.ground.Y_vec;

for i = 1:opts.M    
    sample.x_cell{i} = sample_ABM(params.t_cell{i},sample.D,sample.X_vec(i),sample.T_vec(i));
    % sample.x_cell{i} = opts.ground.x_cell{i};
    sample.y_cell{i} = sample_ABM(params.t_cell{i},sample.D,sample.Y_vec(i),sample.T_vec(i));
    % sample.y_cell{i} = opts.ground.y_cell{i};
end

%% book-keeping
sample.P_vec = get_log_probs(sample.v,sample.D,sample.T_vec, ...
                                    sample.X_vec,sample.Y_vec, ...
                                    sample.x_cell,sample.y_cell,params);

sample.rec = repmat([0;realmin;realmin],1,1);

