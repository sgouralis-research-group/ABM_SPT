
figure(21)
set(gcf,'Name','MCMC')
set(gcf,'windowstyle','docked')
clf

% define burn-in
i_cut = floor(0.25*double(chain.i(end)));
i_lim = [i_cut chain.i(end)+1];
j_idx = find(chain.i >= i_cut);


tiledlayout(1,1)

ax = nexttile;
histogram2(chain.D(j_idx),chain.v(j_idx),'normalization','pdf','ShowEmptyBins','on','DisplayStyle','tile','EdgeColor','none')
xlabel(['D (',chain.params.units.dist,'^2/',chain.params.units.time,')'] )
ylabel(['v (',chain.params.units.dist,'^2)'] )


aD = nexttile('north');
histogram(chain.D(j_idx),'normalization','pdf','FaceAlpha',1,'EdgeColor','none')
xlim(ax.XLim)
box off
aD.YAxis.Visible='off';
aD.XTick = aD.XTick;
aD.XTickLabel = [];
aD.Color = 'none';


av = nexttile('east');
histogram(chain.v(j_idx),'orientation','horizontal','normalization','pdf','FaceAlpha',1,'EdgeColor','none')
ylim(ax.YLim)
box off
av.XAxis.Visible='off';
av.YTick = av.YTick;
av.YTickLabel = [];
av.Color = 'none';


%% add prior
D_lin = linspace(aD.XLim(1),aD.XLim(2),1e2);
line(aD,D_lin,get_inv_gamma(D_lin,chain.params.D_prior_A,(chain.params.D_prior_A-1)*chain.params.D_prior_ref),'color','m','linewidth',2)

v_lin = linspace(av.YLim(1),av.YLim(2),1e2);
line(av,get_inv_gamma(v_lin,chain.params.v_prior_A,(chain.params.v_prior_A-1)*chain.params.v_prior_ref),v_lin,'color','m','linewidth',2)

%% add ground
if isfield( chain.params,'ground' )
    xline(ax,chain.params.ground.D,'-','color','g','linewidth',2)
    yline(ax,chain.params.ground.v,'-','color','g','linewidth',2)
    xline(aD,chain.params.ground.D,'-','color','g','linewidth',2)
    yline(av,chain.params.ground.v,'-','color','g','linewidth',2)
end

%%
linkaxes([ax aD],'x')
linkaxes([ax av],'y')

%% inv-gamma density
function p = get_inv_gamma(t,A,B)
t = max(t,realmin);
p = exp( A*log(B)-gammaln(A)-(A+1)*log(t)-B./t );
end
