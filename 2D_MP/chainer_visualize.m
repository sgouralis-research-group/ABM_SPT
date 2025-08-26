function Gim = chainer_visualize(Gim,chain)
%% init
if isempty(Gim)
    
    num = 4;
    mum = 3;
    
    chain_i = double(chain.i(1)) + chain.stride*(0:chain.length-1)';
    i_lim = [max(chain_i(1),0.01*chain_i(end)) chain_i(end)+1];

    Gim.t_demo_min = chain.params.T_prior_min - 0.25*(chain.params.T_prior_max-chain.params.T_prior_min);
    Gim.t_demo_max = chain.params.T_prior_max + 0.25*(chain.params.T_prior_max-chain.params.T_prior_min);


    figure(10)
    set(gcf,'windowstyle','docked')
    clf

    tiledlayout(num,mum,'TileSpacing','compact','Padding','compact')
     

    % --- Sample ----------------------------------------------------------
    ax_x = nexttile(0*mum+1,[2 mum-1]);
    ax_y = nexttile(2*mum+1,[2 mum-1]);

    line(ax_x,chain.params.t,chain.params.wx,'marker','o','color','k','linestyle','none')
    line(ax_y,chain.params.t,chain.params.wy,'marker','o','color','k','linestyle','none')
    Gim.x_demo = line(ax_x,nan(2,1),nan(2,1),'color','c');
    Gim.y_demo = line(ax_y,nan(2,1),nan(2,1),'color','c');
    Gim.x      = line(ax_x,chain.params.t,nan(chain.params.N,1),'marker','.','linestyle','none');
    Gim.y      = line(ax_y,chain.params.t,nan(chain.params.N,1),'marker','.','linestyle','none');
    Gim.X      = line(ax_x,nan(1),nan(1),'color','m','marker','*','linestyle','none');
    Gim.Y      = line(ax_y,nan(1),nan(1),'color','m','marker','*','linestyle','none');
    xlim(ax_x,[Gim.t_demo_min Gim.t_demo_max])
    xlim(ax_y,[Gim.t_demo_min Gim.t_demo_max])
    ylim(ax_x,chain.params.X_prior_ref + 5*sqrt(chain.params.X_prior_U)*[-1 +1])
    ylim(ax_y,chain.params.Y_prior_ref + 5*sqrt(chain.params.Y_prior_U)*[-1 +1])
    xline(ax_x,chain.params.T_prior_min,'--','label','T_{min}','LabelOrientation','horizontal','LabelHorizontalAlignment','right')
    xline(ax_y,chain.params.T_prior_min,'--','label','T_{min}','LabelOrientation','horizontal','LabelHorizontalAlignment','right')
    xline(ax_x,chain.params.T_prior_max,'--','label','T_{max}','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    xline(ax_y,chain.params.T_prior_max,'--','label','T_{max}','LabelOrientation','horizontal','LabelHorizontalAlignment','left')
    ylabel(ax_x,['space x (',chain.params.units.dist,')'])
    ylabel(ax_y,['space y (',chain.params.units.dist,')'])

    xline(ax_x,chain.params.t,':')
    xline(ax_y,chain.params.t,':')
    yline(ax_x,chain.params.X_prior_ref,'--','Label','X_{ref}','LabelHorizontalAlignment','left')
    yline(ax_y,chain.params.Y_prior_ref,'--','Label','Y_{ref}','LabelHorizontalAlignment','left')
    yline(ax_x,chain.params.X_prior_ref+3*sqrt(chain.params.X_prior_U),':','Label','X_{ref}+3\surd{U_X}','LabelHorizontalAlignment','left')
    yline(ax_y,chain.params.Y_prior_ref+3*sqrt(chain.params.Y_prior_U),':','Label','Y_{ref}+3\surd{U_Y}','LabelHorizontalAlignment','left')
    yline(ax_x,chain.params.X_prior_ref-3*sqrt(chain.params.X_prior_U),':','Label','X_{ref}-3\surd{U_X}','LabelHorizontalAlignment','left')
    yline(ax_y,chain.params.Y_prior_ref-3*sqrt(chain.params.Y_prior_U),':','Label','Y_{ref}-3\surd{U_Y}','LabelHorizontalAlignment','left')

    xlabel(ax_x,['time t (',chain.params.units.time,')'])

    title(ax_x,'MCMC sample')
    subtitle(ax_x,' ')


    % --- MCMC ------------------------------------------------------------
    ax_P = nexttile(1*mum,[2 1]);
    ax_v = nexttile(3*mum,[1 1]);
    ax_D = nexttile(4*mum,[1 1]);
    
    ax_P.YAxisLocation = 'Right';
    ax_v.YAxisLocation = 'Right';
    ax_D.YAxisLocation = 'Right';

    ax_P.XLim = i_lim;
    ax_v.XLim = i_lim;
    ax_D.XLim = i_lim;
    
    title(ax_P,'MCMC chain')
    subtitle(ax_P,['(stride=',num2str(chain.stride),')'])

    xlabel(ax_D, 'MCMC iteration (i)')

    ax_P.YGrid = 'on';
    ax_v.YGrid = 'on';
    ax_D.YGrid = 'on';
  
    ylabel(ax_P,'logP_{post} (nat)' )
    ylabel(ax_v,['v (',chain.params.units.dist,'^2)'] )
    ylabel(ax_D,['D (',chain.params.units.dist,'^2/',chain.params.units.time,')'] )
 
    if ~isempty( chain.ledger )
        xline(ax_P,chain.ledger(end,1));
        xline(ax_v,chain.ledger(end,1));
        xline(ax_D,chain.ledger(end,1));
    end

    Gim.P = line(ax_P,chain_i,chain.P(1,:),'marker','.');
    Gim.v = line(ax_v,chain_i,chain.v     ,'marker','.');
    Gim.D = line(ax_D,chain_i,chain.D     ,'marker','.');


    % --- ground ----------------------------------------------------------
    
    if isfield( chain.params,'ground' )
         line(ax_x,chain.params.t,chain.params.ground.x,'color','g','marker','.');
         line(ax_y,chain.params.t,chain.params.ground.y,'color','g','marker','.');
         line(ax_x,chain.params.ground.T,chain.params.ground.X,'color','g','marker','*','linestyle','none');
         line(ax_y,chain.params.ground.T,chain.params.ground.Y,'color','g','marker','*','linestyle','none');
        yline(ax_v,chain.params.ground.v,'-','color','g')
        yline(ax_D,chain.params.ground.D,'-','color','g')
    end
    
end % init




if chain.params.N<=10
    mult = 10;
elseif chain.params.N<=100
    mult = 2;
else
    mult = 1;
end

 x_demo         = extend_ABM(linspace(Gim.t_demo_min,Gim.t_demo_max,mult*chain.params.N)',chain.sample.x,chain.params.t,chain.sample.D,chain.sample.X,chain.sample.T);
[y_demo,t_demo] = extend_ABM(linspace(Gim.t_demo_min,Gim.t_demo_max,mult*chain.params.N)',chain.sample.y,chain.params.t,chain.sample.D,chain.sample.Y,chain.sample.T);
Gim.x_demo.XData = t_demo;
Gim.y_demo.XData = t_demo;
Gim.x_demo.YData = x_demo;
Gim.y_demo.YData = y_demo;

Gim.x.YData = chain.sample.x;
Gim.y.YData = chain.sample.y;
Gim.X.XData = chain.sample.T;
Gim.Y.XData = chain.sample.T;
Gim.X.YData = chain.sample.X;
Gim.Y.YData = chain.sample.Y;


Gim.P.YData = chain.P(1,:);
Gim.v.YData = chain.v;
Gim.D.YData = chain.D;


drawnow


end % visualize
