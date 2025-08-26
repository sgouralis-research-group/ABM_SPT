

figure(22)
set(gcf,'Name','MCMC')
set(gcf,'windowstyle','docked')
clf

% define burn-in
i_cut = floor(0.25*double(chain.i(end)));
i_lim = [i_cut chain.i(end)+1];
j_idx = find(chain.i >= i_cut);


tiledlayout(4,1,'TileSpacing','compact')

a1 = nexttile;
line(chain.i(j_idx),chain.P_vec(1,j_idx),'marker','.')
ylabel('log P_{post} (nat)')
xlabel(['D (',chain.params.units.dist,'^2/',chain.params.units.time,')'] )

a2 = nexttile([2,1]);
line(chain.i(j_idx),chain.D(j_idx),'marker','.')
ylabel(['D (',chain.params.units.dist,'^2/',chain.params.units.time,')'] )

a3 = nexttile();
line(chain.i(j_idx),chain.v(j_idx),'marker','.')
ylabel(['v (',chain.params.units.dist,'^2)'] )


%% add ground
if isfield( chain.params,'ground' )
    yline(a2,chain.params.ground.D,'-','color','g','linewidth',2)
    yline(a3,chain.params.ground.v,'-','color','g','linewidth',2)
end

%%
linkaxes([a1 a2 a3],'x')




