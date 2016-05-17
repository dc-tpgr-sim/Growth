function [ m_rdf, paravals ] = perp_sum(in_dist, in_map, paravals)
% Takes a 12xN distribution of cells in 2-D perp/para space, and sums over 
% perp values to return a 10xM 1-D distribution.

    % Input 12xM
    % (En, vmag, vpar, vper, pa)_top
    % (En, vmag, vpar, vper, pa)_bottom
    % time, N

    % Input 10xM
    % (En, vmag, vpar, pa)_top
    % (En, vmag, vpar, pa)_bottom
    % time, N

    v_perp = in_dist(in_map.bot.v.perp,:);
    v_para = in_dist(in_map.bot.v.para,:);
    v_N = in_dist(in_map.dist,:);


%    plot(paravals,0:length(paravals)-1,'.')
    n_para = length(paravals);

    [ para_widths, para_deltas ] = half_deltas(paravals);

    m_rdf = zeros(n_para,1);
    parfor i=1:n_para
        para = paravals(i);
        deltas = para_deltas(i:i+1);

        parainds = find(v_para >= para-deltas(1) & v_para < para+deltas(2));
        [ t_perpvals, t_perpinds ] = uniquetol(v_perp(parainds), ...
            0.00001, 'OutputAllIndices', true);
        n_perp = length(t_perpvals);
%        display(num2str(n_perp))

        if n_perp > 1
            % flatten the lists, making a list of all vals, 
            % and a list of indices
            perpvals = [];
            for k=1:n_perp
                perpvals = [ perpvals repmat(t_perpvals(k),1,length(t_perpinds{k})) ];
            end
            perpinds = vertcat(t_perpinds{:});
            
            [ s_perpvals, si_perpvals ] = sort(perpvals);
            
            % indices within this batch of parainds
            si_perpinds = perpinds(si_perpvals);
            % values of f(perp,para)
            s_distN = v_N(parainds(si_perpinds));
            
            % trapezoidal rule function,
            % 1/2 sum( (v_{i+1}-v_i)*(f(v_{i+1})+f(v_i))*v_i )
            delta_v = diff(s_perpvals);
            f_sums = s_distN(1:end-1) + s_distN(2:end);
            f = 0.5*sum( delta_v .* f_sums .* s_perpvals(1:end-1) );

        elseif n_perp == 1
            if length(t_perpinds) > 1
                display('Only one perp value, but multiple indices.  This really shouldn''t happen!')
            end
            f = sum(v_N(parainds(t_perpinds{1})));
            
        else
            f = 0;
            
        end
%        display(num2str(f))

        m_rdf(i) = f;
    end
    
end % parper_rdf()