function [chain] = run_chainer(numTraj, len_t_global, flag_anchor)
    %clear
    format compact
    clc
    
    save_file = [mfilename,'.mat']; % name of the file to store results, by default the same as the filename of this script
    
        %% get data
        [wx_cell,wy_cell,t_cell,t_global,hxx_cell,hyy_cell,hxy_cell, N_vec, M, units_dist,units_time,ground] = generate_synthetic_data(false, numTraj, len_t_global);
        
        %% prepare input for mcmc sampler
        opts.units_time = units_time;
        opts.units_dist = units_dist;
        opts.wx_cell  = wx_cell ;
        opts.wy_cell  = wy_cell ;
        opts.t_cell   = t_cell  ;
        opts.hxx_cell = hxx_cell;
        opts.hyy_cell = hyy_cell;
        opts.hxy_cell = hxy_cell;
        opts.M = M;
        opts.t_global=t_global;
        opts.N_vec = N_vec;
        
        if exist('ground','var')
            opts.ground = ground;
        end
    
        %% init chain
        chain = chainer_main([],0,opts,false,[], flag_anchor);
        
        %clear opts % clear leftovers
        
        chain_temp = chain; % store a temporary copy for debuging
        
        save(save_file,'chain','save_file') % save chain for future processing
        %disp(['SAVED: ', save_file])
        
        
        %% expand chain 
        run_expander
    
    beep
    
    % chain = chainer_main(chain,[],[],true,[],true); % resets chain's record
    
    chain = chainer_main(chain,100,[],true,true); % expands chain's history
    
    % if chain.sizeGB > 1.0
    %     chain = chainer_main(chain,-fix(chain.length/2),[],true,[],[]); % thins chain's history if chain is too big
    % end
    
    save(save_file,'chain','save_file') % save chain for future processing
    disp(['SAVED: ', save_file])

end
