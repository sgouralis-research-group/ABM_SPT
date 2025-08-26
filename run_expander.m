%% expand
while chain.i<1500
    if chain.length < 2
        flag_visual = true;
        d_length =  25;
    else
        flag_visual = false;
        d_length = 250;
    end

    if chain.length >= 1e6
        break
    end

    flag_visual = true;

    % reset records
    chain = chainer_main(chain,      [],[],false,[] , flag_anchor);

    % expand chain
    chain = chainer_main(chain,d_length,[],true,false,flag_anchor);

    % thin chain
    if chain.sizeGB > 1.0
        chain = chainer_main(chain,-fix(chain.length/2),[],true,[], flag_anchor);
    end

    % save results
    save(save_file,'chain','save_file')
    disp(['SAVED: ', save_file])

    % show results
    % if ~flag_visual
    %     Gim = chainer_visualize([],chain);
    % end

    % show_MCMC_results
    % show_MCMC_traces

    % make some noise
    beep
end

