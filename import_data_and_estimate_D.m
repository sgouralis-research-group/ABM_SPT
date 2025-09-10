function [] = import_data_and_estimate_D(filename)
        
data=importdata(filename);

units_time = extractBefore(data.textdata{1},',');
units_dist = extractAfter(data.textdata{1}, ',');

M = max(data.data(:,1));

t_cell   = cell(M,1);
wx_cell  = cell(M,1);
wy_cell  = cell(M,1);
hxx_cell = cell(M,1);
hyy_cell = cell(M,1);
hxy_cell = cell(M,1);
N_vec    = NaN(M,1);

for j = 1:M
    idx = find(data.data(:,1) == j);
    t_cell{j}   = data.data(idx,2);
    N_vec(j)    = size(t_cell{j},1);
    wx_cell{j}  = data.data(idx,3);
    wy_cell{j}  = data.data(idx,4);
    hxx_cell{j} = data.data(idx,5);
    hyy_cell{j} = data.data(idx,6);
    hxy_cell{j} = data.data(idx,7);
end

t_min = min(cellfun(@min,t_cell));
t_max = max(cellfun(@max, t_cell));
t_exp = t_cell{1}(2,:) - t_cell{1}(1,:);

t_global = t_min:t_exp:t_max;

% length of chain
length = 1000;

addpath('2D_MP_v2')
[D_mean,D_std,D_median, D_iqr,chain] = get_D_estimate(M, t_cell, t_global, wx_cell, wy_cell, hxx_cell, ...
                                                       hyy_cell,hxy_cell, units_time, units_dist , true, length);

disp('================================')
disp(['  mean D = ', num2str(D_mean)  , ' ', units_dist,'^2/', units_time])
disp(['   std D = ', num2str(D_std)   , ' ', units_dist,'^2/', units_time])
disp(['median D = ', num2str(D_median), ' ', units_dist,'^2/', units_time])
disp(['   iqr D = ', num2str(D_iqr)   , ' ', units_dist,'^2/', units_time])
disp('================================')
