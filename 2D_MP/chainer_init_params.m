function params = chainer_init_params(opts, flag_anchor)

%% setup
params.units.time = opts.units_time; % units of time
params.units.dist = opts.units_dist; % units of length

params.h_tol = 0.15*min(cellfun(@min, cellfun(@diff, opts.t_cell,'uniformoutput',0)));

params.M = opts.M;
params.t_global=opts.t_global;

for i = 1:params.M
    %% data and dimension checks
    params.N_vec(i) = size(opts.t_cell{i},1);

    params.wx_cell{i}  = reshape(opts.wx_cell{i}, params.N_vec(i), 1);
    params.wy_cell{i}  = reshape(opts.wy_cell{i}, params.N_vec(i), 1);
    params.t_cell{i}   = reshape(opts.t_cell{i}, params.N_vec(i), 1);
    
    % instead save cholesky decomp of h
    params.hxx_cell{i} = reshape(opts.hxx_cell{i}, params.N_vec(i), 1);
    params.hyy_cell{i} = reshape(opts.hyy_cell{i}, params.N_vec(i), 1);
    params.hxy_cell{i} = reshape(opts.hxy_cell{i}, params.N_vec(i), 1);
    

    params.h_n_cell{i} = NaN(2,2,params.N_vec(i));
    params.h_n_cell{i}(1,1,:) = params.hxx_cell{i};
    params.h_n_cell{i}(1,2,:) = params.hxy_cell{i};
    params.h_n_cell{i}(2,1,:) = params.hxy_cell{i};
    params.h_n_cell{i}(2,2,:) = params.hyy_cell{i};

    for j = 1:params.N_vec(i)      
        params.L_chol_cell{i}(:,:,j) = chol(params.h_n_cell{i}(:,:,j), 'lower');
    end

    [params.H_svd_U_cell{i}, params.H_svd_S_cell{i}, ~] ...
              = pagesvd(params.h_n_cell{i},"vector");
    params.h_11_U_cell{i} = squeeze(params.H_svd_U_cell{i}(1,1,:));
    params.h_12_U_cell{i} = squeeze(params.H_svd_U_cell{i}(1,2,:));
    params.h_21_U_cell{i} = squeeze(params.H_svd_U_cell{i}(2,1,:));
    params.h_22_U_cell{i} = squeeze(params.H_svd_U_cell{i}(2,2,:));
    
    params.h_1_S_cell{i} = squeeze(params.H_svd_S_cell{i}(1,1,:));
    params.h_2_S_cell{i} = squeeze(params.H_svd_S_cell{i}(2,1,:));
    
    %% ensure time is ordered
    [~,idx] = sort(params.t_cell{i});

    params.t_cell{i}   = params.t_cell{i}  (idx);
    params.wx_cell{i}  = params.wx_cell{i} (idx);
    params.hxx_cell{i} = params.hxx_cell{i}(idx);
    params.hyy_cell{i} = params.hyy_cell{i}(idx);
    params.hxy_cell{i} = params.hxy_cell{i}(idx);

    
    if flag_anchor 

         %% anchor
         params.T_prior_min_vec(i,1) = min(params.t_cell{i})- 0.1*(max(params.t_cell{i})-min(params.t_cell{i}));
         params.T_prior_max_vec(i,1) = max(params.t_cell{i})+ 0.1*(max(params.t_cell{i})-min(params.t_cell{i}));
        
         params.X_prior_ref_vec(i,1) = mean(params.wx_cell{i});
         params.Y_prior_ref_vec(i,1) = mean(params.wy_cell{i});
         params.X_prior_U_vec(i,1) = 10*((max(params.wx_cell{i})-min(params.wx_cell{i}))^2);   % [dist]^2
         params.Y_prior_U_vec(i,1) = 10*((max(params.wy_cell{i})-min(params.wy_cell{i}))^2);
    else
        %% "no" anchor
        params.T_prior_min_vec(i,1) = min(params.t_cell{i});
        params.T_prior_max_vec(i,1) = min(params.t_cell{i})+1.05*params.h_tol;
        params.X_prior_ref_vec(i,1) = params.wx_cell{i}(1);
        params.Y_prior_ref_vec(i,1) = params.wy_cell{i}(1);
        params.X_prior_U_vec(i,1) = eps;
        params.Y_prior_U_vec(i,1) = eps;
    end
end

params.v_prior_A   = 3;
params.v_prior_ref = mean([mean(cellfun(@var,params.wx_cell)), mean(cellfun(@var,params.wy_cell))])/10; % [dist]^2
    
params.D_prior_A  = 3;
params.D_prior_ref = (mean([mean(cellfun(@var,params.wx_cell)), mean(cellfun(@var,params.wy_cell))])...
                              /mean(cellfun(@range,params.t_cell)));  % [dist]^2/[time]

%% store ground, if exists
if isfield(opts,'ground')
    params.ground = opts.ground;
end


