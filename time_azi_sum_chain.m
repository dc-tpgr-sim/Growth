function [ m_EVPN, smap_EVPN ] = time_azi_sum_chain(in_time, in_dist, in_map, dist_def)

    t_dphi = 3*pi/256; % delta for co-latitude
    t_domega = 0.001; % delta for solid angle in steradians

    % generate dist at top-time t
%    t_dist = dynamic_distribution(in_time, in_dist([in_map.bot.v.mag in_map.bot.alpha],:));
    t_dist = dynamic_distribution(in_time, in_dist, in_map, dist_def);
    m_EVPAN = [ in_dist ; t_dist ]; % Tack distribution on to the rest
    smap_EVPAN = in_map;
    smap_EVPAN.dist = 14;

    % azi_sum
    [ m_EVPN, smap_EVPN ] = azi_sum(m_EVPAN, smap_EVPAN, t_dphi, t_domega);
    
    
end