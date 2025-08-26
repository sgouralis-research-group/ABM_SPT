function D = update_D(x_cell,y_cell,X_vec,Y_vec,T_vec,params)
s=0;
for m=1:params.M
    P = find(T_vec(m) <= params.t_cell{m},1,'first')-1;

    % handle cases where T is outside (t_1, t_N^m)
    if isempty(P)
        P = params.N_vec(m)-1;
        Q = params.N_vec(m);
    elseif P == 0
        P = 1;
        Q = 2;
    else
        Q = P+1;
    end

    s = s + sum( (diff(x_cell{m}(1:P)).^2 ... % pre anchor sum
              +   diff(y_cell{m}(1:P)).^2) ...
              ./ diff(params.t_cell{m}(1:P)) ) ...
           + sum( (diff(x_cell{m}(Q:params.N_vec(m))).^2 ... % post anchor sum
               +   diff(y_cell{m}(Q:params.N_vec(m))).^2) ...
              ./ diff(params.t_cell{m}(Q:params.N_vec(m))) ) ...
           + ((x_cell{m}(P)-X_vec(m))^2 ... % directly left of anchor
           + (y_cell{m}(P)-Y_vec(m))^2) ...
            / (T_vec(m)-params.t_cell{m}(P)) ...
           + ((x_cell{m}(Q)-X_vec(m))^2 ... % directly right of anchor
           + (y_cell{m}(Q)-Y_vec(m))^2) ...
           / (params.t_cell{m}(Q)-T_vec(m));
 end

% posterior sample
D = ((params.D_prior_A-1) * params.D_prior_ref + s/4) ...
    /randg(params.D_prior_A + sum(params.N_vec));

