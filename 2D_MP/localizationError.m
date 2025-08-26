function [locError] = localizationError(chain)
    % calculates localization error of a run

    % remove burn in
    idx=find(chain.i>=0.5*chain.i(end));
    
    chain_dis=NaN(chain.params.M,length(idx));
   
    for m=1:chain.params.M
        chain_dis(m,:)=mean(sqrt((chain.x_cell{m}(:,idx)-chain.params.ground.x_cell{m}).^2 ...
                               + (chain.y_cell{m}(:,idx)-chain.params.ground.y_cell{m}).^2),1);
    end
    
    locError = mean(chain_dis,1);
end
