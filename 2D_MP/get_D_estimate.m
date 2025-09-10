function [D_mean, D_std, D_median, D_iqr, chain] = get_D_estimate(M, t_cell, t_global, wx_cell, wy_cell, ...
                                               hxx_cell, hyy_cell, hxy_cell,...
                                               units_time, units_dist, flag_anchor, ...
                                               chain_length)
    % function get_D_estimate that gets an estimate of D from all 
    % localization data
    
    % get units
    opts.units_time = units_time;
    opts.units_dist = units_dist;

    % get data
    opts.M = M;

    opts.t_global = t_global;

    opts.t_cell = t_cell;

    opts.wx_cell = wx_cell;
    opts.wy_cell = wy_cell;

    opts.hxx_cell = hxx_cell;
    opts.hyy_cell = hyy_cell;
    opts.hxy_cell = hxy_cell;

    %% init chain
    chain = chainer_main([],0,opts,true,[], flag_anchor);
    
    %% expand chain 
    chain = chainer_main(chain, chain_length, [], true, true, flag_anchor);
    
    %% remove burn
    idx = ceil(chain.length/2) : chain.length ;
    
    D_mean   = mean(chain.D(idx));
    D_std    = std(chain.D(idx));
    D_median = median(chain.D(idx));
    D_iqr    = iqr(chain.D(idx));
    % make noise when done
    beep
end