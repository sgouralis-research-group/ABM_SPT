function [] = read_table(filename)
        
data=importdata(filename);

M =  max(data.data(:,1));

t_cell   = cell(M,1);
wx_cell  = cell(M,1);
wy_cell  = cell(M,1);
hxx_cell = cell(M,1);
hyy_cell = cell(M,1);
hxy_cell = cell(M,1);
N_vec    = NaN(M,1);

for j = 1:M
    I = data.data(:,1) == j;
    t_cell{j}   = data.data(I,2);
    N_vec(j)    = size(t_cell{j},1);
    wx_cell{j}  = data.data(I,3);
    wy_cell{j}  = data.data(I,4);
    hxx_cell{j} = data.data(I,5);
    hyy_cell{j} = data.data(I,6);
    hxy_cell{j} = data.data(I,7);
end

t_min = min(cellfun(@min,t_cell));
t_max = max(cellfun(@max, t_cell));
t_exp = t_cell{1}(2,:) - t_cell{1}(1,:);

t_global = t_min:t_exp:t_max;

% length of chain
length = 1000;

addpath('2D_MP_v2')
[D,var,chain] = get_D_estimate(M, t_cell, t_global, wx_cell, wy_cell, hxx_cell, ...
                           hyy_cell,hxy_cell, 'sec', 'micron' , true, length)
