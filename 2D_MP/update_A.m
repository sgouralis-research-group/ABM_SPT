function [sample_T, X, Y,rec] = update_A( ...
          sample_T, x, y, D, m, params,rec)
    %% sample T
    % slice sampling
    [sample_log_p, mu_prime_x, ups_prime_x, mu_prime_y, ups_prime_y ] = get_log_target(sample_T, D, x, y, m, params);
    
    rec(3)=rec(3)+1;
    for i = 1:poissrnd(5)
        log_uprop = log(rand);
        % define interval 
        Tmin = params.T_prior_min_vec(m);
        Tmax = params.T_prior_max_vec(m);
        while true
            % uniform sample from interval
            proposal = Tmin + (Tmax-Tmin)*rand;
            if min(abs(proposal - params.t_cell{m})) < params.h_tol
                prop_log_p = -inf;
            else
                rec(2) = rec(2)+1;
                % calculate log probability of sample
                [prop_log_p, mu_prime_x, ups_prime_x, mu_prime_y, ups_prime_y] ...
                                 = get_log_target(proposal, D, x, y, m, params);
            end
    
            rec(3)=rec(3)+1;
    
            if ~get_sanity_check(prop_log_p-sample_log_p)
                keyboard
            end
    
            if log_uprop < prop_log_p-sample_log_p
                % Accept
                sample_T = proposal;
                sample_log_p = prop_log_p;
                rec(1)=rec(1)+1;
                break
            else
                % update interval
                if proposal < sample_T
                    Tmin = proposal;
                else
                    Tmax = proposal;
                end
            end
        end
    end

    %% recover X and Y
    X = mu_prime_x + sqrt(ups_prime_x) * randn;
    Y = mu_prime_y + sqrt(ups_prime_y) * randn;
end

function [log_P, mu_prime_x, ups_prime_x, mu_prime_y, ups_prime_y ] = get_log_target(T, D, x, y, m, params)
 
    I = find(T <= params.t_cell{m}, 1, 'first');

    % T < t_1
    if I == 1
        ups_prime_x = 1 / (1/(2*D*(params.t_cell{m}(1)-T)) + 1/params.X_prior_U_vec(m));
        mu_prime_x = (x(1)/(2*D*(params.t_cell{m}(1)-T)) + params.X_prior_ref_vec(m)/params.X_prior_U_vec(m)) ...
                    * ups_prime_x;

        ups_prime_y = 1 / (1/(2*D*(params.t_cell{m}(1)-T))+ 1/params.Y_prior_U_vec(m));
        mu_prime_y = (y(1)/(2*D*(params.t_cell{m}(1)-T)) + params.Y_prior_ref_vec(m)/params.Y_prior_U_vec(m)) ...
                    * ups_prime_y;

        log_J_x = 0.5 * (log(ups_prime_x) - log(params.X_prior_U_vec(m)) - log(2*D*(params.t_cell{m}(1)-T)) ...
                        + (mu_prime_x - x(1))*x(1) / (2*D*(params.t_cell{m}(1)-T)) ...
                        + (mu_prime_x-params.X_prior_ref_vec(m))*params.X_prior_ref_vec(m)/params.X_prior_U_vec(m));
        log_J_y = 0.5 * (log(ups_prime_y) - log(params.Y_prior_U_vec(m)) - log(2*D*(params.t_cell{m}(1)-T)) ...
                        + (mu_prime_y - y(1))*y(1) / (2*D*(params.t_cell{m}(1)-T)) ...
                        + (mu_prime_y-params.Y_prior_ref_vec(m))*params.Y_prior_ref_vec(m)/params.Y_prior_U_vec(m));

     % t_n < T
     elseif isempty(I)
         ups_prime_x = 1 / (1/(2*D*(T-params.t_cell{m}(end)))+ 1/params.X_prior_U_vec(m));
         mu_prime_x = (x(end)/(2*D*(T-params.t_cell{m}(end))) + params.X_prior_ref_vec(m)/params.X_prior_U_vec(m)) ...
                    * ups_prime_x;

        ups_prime_y = 1 / (1/(2*D*(T-params.t_cell{m}(end)))+ 1/params.Y_prior_U_vec(m));
        mu_prime_y = (y(end)/(2*D*(T-params.t_cell{m}(end))) + params.Y_prior_ref_vec(m)/params.Y_prior_U_vec(m)) ...
                    * ups_prime_y;

        log_J_x = 0.5 * (log(ups_prime_x) - log(params.X_prior_U_vec(m)) - log(2*D*(T-params.t_cell{m}(end))) ...
                        + (mu_prime_x - x(end))*x(end) / (2*D*(T-params.t_cell{m}(end))) ...
                        + (mu_prime_x-params.X_prior_ref_vec(m))*params.X_prior_ref_vec(m)/params.X_prior_U_vec(m));
        log_J_y = 0.5 * (log(ups_prime_y) - log(params.Y_prior_U_vec(m)) - log(2*D*(T-params.t_cell{m}(end))) ...
                        + (mu_prime_y - y(end))*y(end) / (2*D*(T-params.t_cell{m}(end))) ...
                        + (mu_prime_y-params.Y_prior_ref_vec(m))*params.Y_prior_ref_vec(m)/params.Y_prior_U_vec(m));
    
    % t_P < T < t_Q
    else
        P = I-1;
        Q = I;
        ups_prime_x = 1 / (1/(2*D*(T-params.t_cell{m}(P))) ...
                          +1/params.X_prior_U_vec(m) ...
                          +1/(2*D*(params.t_cell{m}(Q)-T)));
        mu_prime_x  = (x(P)/(2*D*(T-params.t_cell{m}(P)))...
                      + params.X_prior_ref_vec(m)/params.X_prior_U_vec(m)...
                      + x(Q)/(2*D*(params.t_cell{m}(Q)-T))) ...
                      * ups_prime_x;

        ups_prime_y = 1 / (1/(2*D*(T-params.t_cell{m}(P))) ...
                          +1/params.Y_prior_U_vec(m) ...
                          +1/(2*D*(params.t_cell{m}(Q)-T)));
        mu_prime_y  = (y(P)/(2*D*(T-params.t_cell{m}(P)))...
                      + params.Y_prior_ref_vec(m)/params.Y_prior_U_vec(m)...
                      + y(Q)/(2*D*(params.t_cell{m}(Q)-T))) ...
                      * ups_prime_y;

        log_J_x = 0.5 *   (log(4*D*(params.t_cell{m}(Q)-params.t_cell{m}(P))) ...
                        + (x(Q)-x(P))^2 / (2*D*(params.t_cell{m}(Q)-params.t_cell{m}(P))) ...
                        + log(ups_prime_x)...
                        - log(2*D*(T-params.t_cell{m}(P))) ...
                        - log(params.X_prior_U_vec(m)) ... 
                        - log(2*D*(params.t_cell{m}(Q)-T)) ...
                        + (mu_prime_x-x(P))*x(P)/(2*D*(T-params.t_cell{m}(P))) ...
                        + (mu_prime_x-params.X_prior_ref_vec(m))*params.X_prior_ref_vec(m)/params.X_prior_U_vec(m) ...
                        + (mu_prime_x-x(Q))*x(Q)/(2*D*(params.t_cell{m}(Q)-T)) );

        log_J_y = 0.5 *   (log(4*D*(params.t_cell{m}(Q)-params.t_cell{m}(P))) ...
                        + (y(Q)-y(P))^2 / (2*D*(params.t_cell{m}(Q)-params.t_cell{m}(P))) ...
                        + log(ups_prime_y)...
                        - log(2*D*(T-params.t_cell{m}(P))) ...
                        - log(params.Y_prior_U_vec(m)) ... 
                        - log(2*D*(params.t_cell{m}(Q)-T)) ...
                        + (mu_prime_y-y(P))*y(P)/(2*D*(T-params.t_cell{m}(P))) ...
                        + (mu_prime_y-params.Y_prior_ref_vec(m))*params.Y_prior_ref_vec(m)/params.Y_prior_U_vec(m) ...
                        + (mu_prime_y-y(Q))*y(Q)/(2*D*(params.t_cell{m}(Q)-T)) );
    end

    log_P = log_J_x + log_J_y;
end


