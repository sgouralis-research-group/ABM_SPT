function sample = sampler_update( sample, params )

%% update counter
sample.i = sample.i + 1;
    % random sweep
    for tag = randperm(4)
        switch tag
            case 1
                for m = 1:params.M
                   [sample.T_vec(m), sample.X_vec(m),  sample.Y_vec(m),                       sample.rec(:,tag)] = update_A( ...
                    sample.T_vec(m), sample.x_cell{m}, sample.y_cell{m}, sample.D, m, params, sample.rec(:,tag));
                end
            case 2
                for m = 1:params.M
                    [sample.x_cell{m}, sample.y_cell{m}]...
                     = update_R(sample.T_vec(m), sample.X_vec(m), sample.Y_vec(m),sample.D, sample.v, m, params);
                end 
            case 3
                 sample.D = update_D(sample.x_cell,sample.y_cell,sample.X_vec,sample.Y_vec,sample.T_vec,params);
            case 4
                 sample.v = update_v(sample.x_cell,sample.y_cell,params);
            otherwise
                error('Unknown sampler requested')
        end % switch
    end % tag


%% book-keeping
sample.P_vec = get_log_probs(sample.v,sample.D,sample.T_vec,sample.X_vec, ...
                             sample.Y_vec,sample.x_cell,sample.y_cell, params);
        