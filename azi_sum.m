function [ m_fperppara, out_map, m_intersect, s_uniques ] = azi_sum(in_dist, in_map, t_dalpha, t_domega)
% Takes a 14xN list of cells in a distribution in 3-D perp/para/azi space,
% and sums over Azimuthseseses to return a 12xM list of 2-D reduced 
% distribution functions in v_perp-v_para space.  Optionally returns its 
% intersection list and a structure of the unique perp and para values.

    % Input 14xN matrix
    % (En, vmag, vpar, vper, pa, azi)_top
    % (En, vmag, vpar, vper, pa, azi)_bottom
    % time, distN

    % Returns 12xM
    % (En, vmag, vpar, vper, pa)_top
    % (En, vmag, vpar, vper, pa)_bottom
    % time, N

    v_perp = in_dist(in_map.bot.v.perp,:);
    v_para = in_dist(in_map.bot.v.para,:);

    % get cell of tuple-matches
    [ m_intersect, s_uniques ] = azi_sum_stash(v_perp, v_para);
    n_cells = length(m_intersect);

    m_fperppara = zeros(12,n_cells);
    for i=1:n_cells
        v_indices = m_intersect{i};

        v_distN = in_dist(in_map.dist,v_indices);        

        % Since these are limited to a single vperp,vpara, they all have
        % the same alpha, i.e. they're in an azimuthal ring.  Because
        % that's exactly how azimuths were defined.  Thus, dtheta is just
        % 2pi/(# of points).  We can just factor that 
        t_dtheta = 2*pi/length(v_indices);        
        t_sumN = sum(v_distN*t_dtheta);

        m_fperppara(:,i) = [ in_dist([ ...
            in_map.top.En in_map.top.v.mag ... 
            in_map.top.v.perp in_map.top.v.para in_map.top.alpha ...
            in_map.bot.En in_map.bot.v.mag ... 
            in_map.bot.v.perp in_map.bot.v.para in_map.bot.alpha ...
            in_map.time],v_indices(1)) ; t_sumN ];
            % The values from the input should be identical for all v_indices()
            
    end

    % create new output field map
    out_map.top.En = 1; out_map.top.v.mag = 2;
    out_map.top.v.perp = 3; out_map.top.v.para = 4;
    out_map.top.alpha = 5;
    out_map.bot.En = 6; out_map.bot.v.mag = 7;
    out_map.bot.v.perp = 8; out_map.bot.v.para = 9;
    out_map.bot.alpha = 10;
    out_map.time = 11; out_map.dist = 12;
    
end