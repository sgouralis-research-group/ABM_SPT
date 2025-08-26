function [wx_cell,wy_cell,t_cell,t_global,hxx_cell,hyy_cell,hxy_cell, N_vec, M ,units_dist,units_time,ground]...
           = generate_synthetic_data(demo, numTraj, len_t_global)

if nargin == 0
    demo = true;
end

% rng(3);

% set units
units_dist = '\mum';
units_time = 's';

% number of particles
M = numTraj;

% set parameters
v = 0.025; % [dist]^2
D = 0.1; % [dist]^2/[time]
ground.v = v;
ground.D = D;

T_min =  2; % [time]
T_max = 18; % [time]

% set assesment points and error bars
t_min = T_min - 0.25*(T_max-T_min);
t_max = T_max + 0.25*(T_max-T_min);

X_ref = +5;  % [dist]
Y_ref = -1;  % [dist]
UX    =  10;  % [dist]^2
UY    =  10;  % [dist]^2

% length of M trajectories
N_vec = NaN(M,1);

% allocate
t_cell = cell(M,1);
T_vec = NaN(M,1);
X_vec = NaN(M,1);
Y_vec = NaN(M,1);
x_cell = cell(M,1);
y_cell = cell(M,1);
wx_cell = cell(M,1);
wy_cell = cell(M,1);
hxx_cell = cell(M,1);
hyy_cell = cell(M,1);
hxy_cell = cell(M,1);

t_global = linspace(t_min, t_max, len_t_global)';

for i = 1:M
    % length_param is geornd param that controls length of treajectories
    % 0.25 -> average length of around 4
    length_param = 0.25;

    n1 = randi(size(t_global,1)-2);
    n2 = t_global(end)+1;
    while n2>=length(t_global)
        n2 = n1+geornd(length_param)+1;
    end

    N1 = min(n1,n2);
    N2 = max(n1,n2);

    N_vec(i) = (N2-N1)+1;

    % % for some reason, the geometric rv version doesnt work when doing
    % % parallel runs. uncomment below for use with parallel chains, comment
    % % above
    % n1 = randi(length(t_global)-1);
    % n2 = randi(length(t_global));
    % 
    % n2 = n2+(n2==n1);
    % N1 = min(n1,n2);
    % N2 = max(n1,n2);
    % 
    % N_vec(i) = (N2-N1)+1;

    hxx_cell{i} =  ones(N_vec(i),1);
    hyy_cell{i} =  ones(N_vec(i),1);
    hxy_cell{i} = zeros(N_vec(i),1);
    
    t_cell{i} = t_global(N1:N2);

    %% generate the anchor and the trajectory
    T_vec(i) = t_cell{i}(1) + (range(t_cell{i}))*rand;
    X_vec(i) = X_ref + sqrt(UX)*randn;
    Y_vec(i) = Y_ref + sqrt(UY)*randn;
    
    x_cell{i} = sample_ABM(t_cell{i},D,X_vec(i),T_vec(i));
    y_cell{i} = sample_ABM(t_cell{i},D,Y_vec(i),T_vec(i));
    
    
    %% add noise
    % Cholesky analytically
    L11 = sqrt(hxx_cell{i});
    L12 = hxy_cell{i}./L11;
    L22 = sqrt(hyy_cell{i}-L12.^2);
    
    d = sqrt(v)*randn(N_vec(i),2);
    wx_cell{i} = x_cell{i} + L11.*d(:,1);
    wy_cell{i} = y_cell{i} + L12.*d(:,1) + L22.*d(:,2);
end

%% maintain ground
ground.T_vec = T_vec;
ground.X_vec = X_vec;
ground.Y_vec = Y_vec;
ground.x_cell = x_cell;
ground.y_cell = y_cell;

%% demo
% pick a trajectory
i = randi(M);
if demo

    % extend the trajectory for clear demonstration
    t_min_demo = t_min - 0.3*(t_max-t_min);
    t_max_demo = t_max + 0.3*(t_max-t_min);
     x_demo         = extend_ABM(linspace(t_min_demo,t_max_demo,100*N_vec(i))',x_cell{i},t_cell{i},D,X_vec(i),T_vec(i));
    [y_demo,t_demo] = extend_ABM(linspace(t_min_demo,t_max_demo,100*N_vec(i))',y_cell{i},t_cell{i},D,Y_vec(i),T_vec(i));


    %% Visualize generated data if desired, uncomment below
    figure(1)
    %set(gcf,'windowstyle','docked')
    clf

    tiledlayout(1,1,'TileSpacing','compact')

    s = sprintf('Simulated ABM of Length %d', N_vec(i) );
    sgtitle(s, 'FontSize', 20,'interpreter','latex')
    nexttile()

    line(T_vec(i),X_vec(i),'color','r','marker','*','linestyle','none')  % shows the anchor
    line(t_demo,x_demo)                                    % shows the trajectory
    line(t_cell{i},x_cell{i},'color','k','marker','.','linestyle','none', 'markersize', 10)  % shows the assesments
    %line(t_cell{i},wx_cell{i},'color','k','marker','o','linestyle','none')  % shows the measurments

    xline(t_cell{i},':')

    legend('$(X, T)$','$x(\cdot)$','$x(t_n)$','$t_n$', ...%'$w_{x,n}$','$t_n$',...
           'Location','northeast','box','on','orientation','vertical','AutoUpdate','Off', ...
           'interpreter','latex')

    %xlim( [t_min t_max] + 0.075*(t_max-t_min)*[-1 +1] )
    %ylim(X_ref + 4*sqrt(UX)*[-1 +1])

    xlabel(['time t (',units_time,')'],'interpreter','latex')
    ylabel('space x ($\mu$m)','interpreter','latex')

    yline(X_vec(i),'--','color', 'r', 'Label','$X$','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','interpreter','latex')
    xline(T_vec(i),'--','color', 'r','Label','$T$','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','interpreter','latex')

    %yline(X_ref,'--','color','b','Label','$X_{ref}$','LabelHorizontalAlignment','left','LabelOrientation','horizontal','interpreter','latex')
    %yline(X_ref+3*sqrt(UX),':','color','b','Label','$X_{ref}+3\sqrt{U_X}$','LabelHorizontalAlignment','left','interpreter','latex')
    %yline(X_ref-3*sqrt(UX),':','color','b','Label','$X_{ref}-3\sqrt{U_X}$','LabelHorizontalAlignment','left','interpreter','latex')

    %xline(t_min,'color', 'b','Label','$t_{min}$','LabelHorizontalAlignment','left','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','interpreter','latex')
    %xline(t_max,'color', 'b','Label','$t_{max}$','LabelHorizontalAlignment','right','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','interpreter','latex')

    %xline(T_min,'color', 'b','Label','$T_{min}$','LabelHorizontalAlignment','left','LabelOrientation','horizontal','interpreter','latex')
    %xline(T_max,'color', 'b','Label','$T_{max}$','LabelHorizontalAlignment','right','LabelOrientation','horizontal','interpreter','latex')


    % nexttile()
    % 
    % line(T_vec(i),Y_vec(i),'color','r','marker','*','linestyle','none')  % shows the anchor
    % line(t_demo,y_demo)                                    % shows the trajectory
    % line(t_cell{i},y_cell{i},'color','k','marker','.','linestyle','none')  % shows the assesments
    % line(t_cell{i},wy_cell{i},'color','k','marker','o','linestyle','none')  % shows the measurments
    % 
    % xline(t_cell{i},':')
    % 
    % legend('$Y$','$y(\cdot)$','$y(t_n)$','$w_{y,n}$','$t_n$',...
    %        'Location','northoutside','box','off','orientation','horizontal','AutoUpdate','Off', ...
    %        'interpreter','latex')
    % 
    % xlim( [t_min t_max] + 0.15*(t_max-t_min)*[-1 +1] )
    % ylim(Y_ref + 5*sqrt(UY)*[-1 +1])
    % 
    % xlabel(['time t (',units_time,')'],'interpreter','latex')
    % ylabel('space y ($\mu$m)','interpreter','latex')
    % 
    % yline(Y_vec(i),'--','color','r','Label','$Y$','LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','interpreter','latex')
    % xline(T_vec(i),'--','color','r','Label','$T$','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','LabelOrientation','horizontal','interpreter','latex')
    % 
    % yline(Y_ref,'--','color','b','Label','$Y_{ref}$','LabelHorizontalAlignment','left','LabelOrientation','horizontal','interpreter','latex')
    % yline(Y_ref+3*sqrt(UY),':','color','b','Label','$Y_{ref}+3\sqrt{U_Y}$','LabelHorizontalAlignment','left','interpreter','latex')
    % yline(Y_ref-3*sqrt(UY),':','color','b','Label','$Y_{ref}-3\sqrt{U_Y}$','LabelHorizontalAlignment','left','interpreter','latex')
    % 
    % xline(t_min,'color','b','Label','$t_{min}$','LabelHorizontalAlignment','left','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','interpreter','latex')
    % xline(t_max,'color','b','Label','$t_{max}$','LabelHorizontalAlignment','right','LabelOrientation','horizontal','LabelVerticalAlignment','bottom','interpreter','latex')
    % 
    % xline(T_min,'color','b','Label','$T_{min}$','LabelHorizontalAlignment','left','LabelOrientation','horizontal','interpreter','latex')
    % xline(T_max,'color','b','Label','$T_{max}$','LabelHorizontalAlignment','right','LabelOrientation','horizontal','interpreter','latex')
    % 
    % fontsize(gcf,scale=1.3)
end % demo

% clear leftovers
if nargout == 0
    clear 
end


