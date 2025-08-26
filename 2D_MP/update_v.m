function v = update_v(x_cell,y_cell,params)
s=0;

for m=1:params.M

    s = sum(pagemldivide(params.L_chol_cell{m}, ...
                 reshape([params.wx_cell{m}, params.wy_cell{m}]-[x_cell{m}, y_cell{m}], ...
                          2,1,[])).^2,'all') ...
        + s;
    
end

% posterior sample
v = (s/2 + (params.v_prior_A-1)*params.v_prior_ref) ...
    / randg(params.v_prior_A + sum(params.N_vec));

