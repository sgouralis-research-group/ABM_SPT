% This code is used to run run_chainer in parallel when you want to do many
% runs to compare performance with trajectory lengths, with/without anchor,
% etc.

% show_MCMC_traces and show_MCMC_results in run_chainer must be 
% commented in order to this script to run, as well as geometric rv generation
% of the trajectory lengths in chain_init_params
% also commenting out the printed
% progress text in run_chainer imporves performance.

numruns = 20;

% Ms = [1 5 15 30  40  50 70 80 100 120];
Ms = round(linspace(5,  120, numruns));

% t_globs = [1000 400 90 50 40 25 22 20 15 10];
% Ms =    [5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, 155, 175, 185, 195];
t_globs = [10 15  20  25  30  35  40  40  40  40  40   40   50   60   70   80  80  90 100 120];

chains = cell(numruns, 1);

bin = [true; false];
for b = 1:2
    parfor i=1:numruns
        % chains{i}=run_chainer(Ms(i), t_globs(i), true);
        chains{i}=run_chainer(Ms(i), t_globs(i), bin(b));
        % disp(i);
    end
    % disp(b);
    % c=[chains{:}];
    c{b} = [chains{:}];
end


