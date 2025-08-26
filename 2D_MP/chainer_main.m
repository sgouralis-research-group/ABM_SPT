function chain = chainer_main(chain_init,d_length,opts,flag_status,flag_visual, flag_anchor)
% to init:
% chain = chainer_main([]   ,  0, opts, true, []  );
% to expand:
% chain = chainer_main(chain,+25, []  , true, true);
% to reduce:
% chain = chainer_main(chain,-10, []  , true, []  );
% to reset record:
% chain = chainer_main(chain, [], []  , true, []  );

tic_id = tic;

% initialize the seed or use seed for last expansion
if isempty(chain_init)
    rng('shuffle');
else
    rng(chain_init.random);
end

% init chain --------------------------------------------------------------
if d_length == 0

    % MCMC
    chain.params = chainer_init_params(opts,flag_anchor);
    chain.length = 1;
    chain.stride = 1;
    chain.ledger = nan(0,2); % total wall time without Initialization
    chain.sizeGB = nan;      % current chain memory size
    chain.record = [];       % acceptance rates
    chain.sample = [];

    chain.sample = chainer_init_sample(chain.params,opts);
    
    % history
    chain.i = cast( chain.sample.i, 'uint64');
    chain.P_vec = cast( chain.sample.P_vec, 'double');
    chain.v = cast( chain.sample.v, 'single');
    chain.D = cast( chain.sample.D, 'single');
    chain.T_vec = cast( chain.sample.T_vec, 'single');
    chain.X_vec = cast( chain.sample.X_vec, 'single');
    chain.Y_vec = cast( chain.sample.Y_vec, 'single');
    chain.x_cell = chain.sample.x_cell;
    chain.y_cell = chain.sample.y_cell;
    
    if flag_status
        disp('CHAINER: chain initiated')
    end

% expand chain ------------------------------------------------------------
elseif d_length > 0
    chain.params = chain_init.params;
    chain.length = chain_init.length + d_length;
    chain.stride = chain_init.stride;
    chain.ledger = chain_init.ledger;
    chain.sizeGB = nan;
    chain.record = chain_init.record;
    chain.sample = chain_init.sample;

    chain.i = cat( 2 , chain_init.i, zeros( 1              , d_length , 'like',chain_init.i ) );
    chain.P_vec = cat( 2 , chain_init.P_vec,   nan( get_log_probs  , d_length , 'like',chain_init.P_vec ) );
    chain.v = cat( 2 , chain_init.v,   nan( 1              , d_length , 'like',chain_init.v ) );
    chain.D = cat( 2 , chain_init.D,   nan( 1              , d_length , 'like',chain_init.D ) );
    chain.T_vec = cat( 2 , chain_init.T_vec,   nan( chain.params.M, d_length , 'like',chain_init.T_vec ) );
    chain.X_vec = cat( 2 , chain_init.X_vec,   nan( chain.params.M, d_length , 'like',chain_init.X_vec ) );
    chain.Y_vec = cat( 2 , chain_init.Y_vec,   nan( chain.params.M, d_length , 'like',chain_init.Y_vec ) );
    
    % each cell i will contain the history of its trajectory
    for i = 1:chain.params.M
        chain.x_cell{i} = cat( 2 , chain_init.x_cell{i},   nan( chain.params.N_vec(i) , d_length , 'like', chain_init.x_cell{i} ) );
        chain.y_cell{i} = cat( 2 , chain_init.y_cell{i},   nan( chain.params.N_vec(i) , d_length , 'like', chain_init.y_cell{i} ) );
    end
    
    % if flag_visual
    %     Gim = chainer_visualize([],chain);
    % end
    
    %---------------------------- expand chain
    r = chain_init.length+1;
    while r <= chain.length
        
        chain.sample = sampler_update(chain.sample,chain.params);

        if mod(chain.sample.i,chain.stride) == 0
            
            chain.i(r) = chain.sample.i;
            chain.P_vec(:,r) = chain.sample.P_vec;
            chain.v(r) = chain.sample.v;
            chain.D(r) = chain.sample.D;
            chain.T_vec(:,r) = chain.sample.T_vec;
            chain.X_vec(:,r) = chain.sample.X_vec;
            chain.Y_vec(:,r) = chain.sample.Y_vec;
            
            for i = 1:chain.params.M
                chain.x_cell{i}(:,r) = chain.sample.x_cell{i};
                chain.y_cell{i}(:,r) = chain.sample.y_cell{i};
            end

            % if flag_visual
            %     chainer_visualize(Gim,chain);
            % end
            
            if flag_status
                disp([  'i = ', num2str(chain.sample.i,'%d'), ...
                     ' - acc = ', ...
                                num2str( chain.sample.rec(1,:)./chain.sample.rec(2,:)  * 100 ,'%#6.2f') , ' %', ...
                     ])
            end
            
            r = r+1;
        end
    end    

    if flag_status
        disp('CHAINER: chain expanded')
    end

% reduce chain ------------------------------------------------------------
elseif d_length < 0

    d_length = min(-d_length,chain_init.length);

    chain.params = chain_init.params;
    chain.length = d_length;
    chain.stride = nan;
    chain.ledger = chain_init.ledger;
    chain.sizeGB = nan;
    chain.record = chain_init.record;
    chain.sample = chain_init.sample;
    
    ind = mod(chain_init.length,d_length)+(floor(chain_init.length/d_length)*(1:d_length));

    chain.i = chain_init.i(ind);
    chain.P_vec = chain_init.P_vec(:, ind);
    chain.v = chain_init.v(ind);
    chain.D = chain_init.D(ind);
    chain.T_vec = chain_init.T_vec(:, ind);
    chain.X_vec = chain_init.X_vec(:, ind);
    chain.Y_vec = chain_init.Y_vec(:, ind);

    for i=1:chain.params.M
        chain.x_cell{i} = chain_init.x_cell{i}(:,ind);
        chain.y_cell{i} = chain_init.y_cell{i}(:,ind);
    end
    
    chain.stride = double(chain.i(2)-chain.i(1));
    
    if flag_status
        disp('CHAINER: chain reduced')
    end
    
    
% reset chain -------------------------------------------------------------
elseif isempty(d_length)

    chain = chain_init;
    
    chain.record = [chain.record; [chain.sample.i,chain.sample.rec(1,:)./chain.sample.rec(2,:)] ];

    chain.sample.rec(1,:) = 0;
    chain.sample.rec(2,:) = realmin;
    chain.sample.rec(3,:) = realmin;

    if flag_status
        disp('CHAINER: chain reset')
    end
    
end


% store the seed for future expansion
chain.random = rng();



%% book-keeping
chain.sizeGB = get_sizeGB(chain);               % mem size

% ledger
wall_time = toc(tic_id);
chain.ledger = [chain.ledger; double(chain.i(end)), wall_time];

if flag_status
    disp(['( wall time = ',num2str(wall_time),' s, overall wall time = ',num2str(sum(chain.ledger(:,2))),' s )'])
end
end

%% auxiliary functions
function sizeGB = get_sizeGB(chain)
    sizeGB = whos( inputname(1) );
    sizeGB = sizeGB.bytes/1024^3;
end

