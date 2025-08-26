function [x,y] = update_R(T,X,Y,D,v,m,params)
    % allocate
    mubars =  NaN(2,    params.N_vec(m));
    nubars =  NaN(2, 2, params.N_vec(m));
    muhats =  NaN(2,    params.N_vec(m));
    nuhats =  NaN(2, 2, params.N_vec(m));
    z =       NaN(2,    params.N_vec(m));

    %% first do t_n > T
    n = find(T<params.t_cell{m}, 1, 'first');
    if isempty(n)
        n = params.N_vec(m);
    end

    mubars(:,n) = [X ; Y];
    nubars(:,:,n) = 2*D * [abs(T-params.t_cell{m}(n)), 0  ; ...
                            0, abs(T-params.t_cell{m}(n))] ;

    G = nubars(:,:,n)/(v*params.h_n_cell{m}(:,:,n)+nubars(:,:,n));
    %svd(G); % singular check

    muhats(:,n) = mubars(:,n) + G*([params.wx_cell{m}(n); params.wy_cell{m}(n)] - mubars(:,n));
    nuhats(:,:,n) = (eye(2)-G)*nubars(:,:,n);
    
    % calculate the filter parameters
    for i = n+1:params.N_vec(m)
        mubars(:,i) = muhats(:,i-1);
        nubars(:,:,i)= step_covariance(i, params.t_cell{m}, D) + nuhats(:,:,i-1);
        
        G = nubars(:,:,i)/(v*params.h_n_cell{m}(:,:,i)+nubars(:,:,i));
        
        muhats(:,i) = mubars(:,i) ...
                    + G ...
                    * ([params.wx_cell{m}(i); params.wy_cell{m}(i)]-mubars(:,i));
    
        nuhats(:,:,i) = (eye(2)-G)*nubars(:,:,i);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % backwards sampling %
    %%%%%%%%%%%%%%%%%%%%%%

    % first sample
    [U,S,~]=svd(nuhats(:,:,params.N_vec(m)));
    z(:,params.N_vec(m)) = muhats(:,params.N_vec(m)) + U*sqrt(S)*randn(2,1);

    for i = params.N_vec(m)-1:-1:n
        J = nuhats(:,:,i) / nubars(:,:,i+1);

        [U,S,~]=svd(nuhats(:,:,i) - J*nubars(:,:,i+1)*J');

        z(:,i) = (muhats(:,i) + J*(z(:,i+1)-mubars(:,i+1))) + U*sqrt(S)*randn(2,1);
    end

    %% Now do t_n-1 < T
    if n == 1
        n = 2;
    end

    mubars(:,n-1) = [X ; Y];
    nubars(:,:,n-1) = 2*D* [abs(T-params.t_cell{m}(n-1)), 0  ; ...
                             0, abs(T-params.t_cell{m}(n-1))] ;
    
    G = nubars(:,:,n-1)/(v*params.h_n_cell{m}(:,:,n-1)+nubars(:,:,n-1));
    %svd(G); % singular check

    muhats(:,n-1) = mubars(:,n-1) + G*([params.wx_cell{m}(n-1); params.wy_cell{m}(n-1)] - mubars(:,n-1));
    nuhats(:,:,n-1) = (eye(2)-G)*nubars(:,:,n-1);

    % calculate the filter parameters (moving backwards)
    for i = n-1:-1:1
        mubars(:,i) = muhats(:,i+1);
        nubars(:,:,i) = step_covariance(i+1, params.t_cell{m}, D) + nuhats(:,:,i+1);
        G = nubars(:,:,i)/(v*params.h_n_cell{m}(:,:,i)+nubars(:,:,i));
        
        muhats(:,i) = mubars(:,i) ...
                    + G ...
                    * ([params.wx_cell{m}(i); params.wy_cell{m}(i)]-mubars(:,i));
    
        nuhats(:,:,i) = (eye(2)-G)*nubars(:,:,i);
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % backwards sampling %
    %%%%%%%%%%%%%%%%%%%%%%
    % (moving forwards)
    % first sample
    [U,S,~]=svd(nuhats(:,:,1));
    z(:,1) = muhats(:,1) + U*sqrt(S)*randn(2,1);

    for i = 2:n-1
        J = nuhats(:,:,i) / nubars(:,:,i-1);
        
        [U,S,~]=svd(nuhats(:,:,i) - J*nubars(:,:,i-1)*J');
        
        z(:,i) = (muhats(:,i) + J*(z(:,i-1)-mubars(:,i-1))) + U*sqrt(S)*randn(2,1);
    end
    
    %% x and y
    x = z(1,:)';
    y = z(2,:)';

end

function V = step_covariance(i,ts,D)
    V = 2*D*[abs(ts(i)-ts(i-1)), 0 ; ...
             0, abs(ts(i)-ts(i-1))];
end
