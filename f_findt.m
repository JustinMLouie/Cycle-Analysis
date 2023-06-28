%% Returns an interpolated time for a percentage pressure drop
function t_int = f_findt(t, p_norm, drop)
    % t = capture time interval
    % p_norm = normalized capture pressure profile
    % drop = percent pressure drop to be interpolated

    % If/else statement to calculate times of multiple drop targets

    if length(drop) <= 1
        % Find difference between normalized profile and target
        target = 1 - drop;
        diff = p_norm - target;

        % Discretize interval of sign flip
        ind_neg = find(diff < 0, 1); % first index of sign flip
        ind_pts = [ind_neg - 1, ind_neg];
        tq = linspace(t(ind_neg - 1), t(ind_neg), 5e2); % new time interval to query
        
        % Repeat process within new interval, take sign flip as solved t
        p_int = interp1(t(ind_pts), p_norm(ind_pts), tq); % interpolated pressures
        t_sign = tq(p_int - target < 0);
        t_int = t_sign(1);
        
    else
        % initialize time array
        t_int = zeros(size(drop));

        % loop through array of pressure drops
        for i = 1:length(drop)
            % find difference between normalized profile, target
            target = (1 - drop(i));
            diff = p_norm - target;

            % discretize interval of sign flip
            ind_neg = find(diff < 0, 1); % first index of sign flip
            ind_pts = [ind_neg-1, ind_neg];
            tq = linspace(t(ind_neg-1), t(ind_neg), 5e2); % new time interval to query

            % repeat process within new interval, take sign flip as solved t
            p_int = interp1(t(ind_pts), p_norm(ind_pts), tq); % interpolated pressures
            t_sign = tq(p_int - target < 0);
            t_int(i) = t_sign(1);

        end
    end
end