function [ m_intersect, s_uniques ] = azi_sum_stash(v_perp,v_para)
% Finds unique (v_perp,v_para) tuples and returns the indices from the data
% that hit those tuples, i.e. a cell of arrays of data poitns with the same
% (v_perp,v_para), but different azimuths.
%
% Will cache this search for inputs which match checksums, because the
% uniquetol()s and intersections are rather slow.

    % uniquetol() and the set intersection stuff are very time consuming, so
    % we'll cache those results and a checksum.
    persistent perpcs paracs sm_intersect ss_uniques

    % First check if we've got an accurate cache.
    newperpcs = array_checksum(v_perp); % checksum of perp velocities
    newparacs = array_checksum(v_para); % checksum of para velocities

    if ~isequal(perpcs,newperpcs) || ~isequal(paracs,newparacs) || isempty(sm_intersect)
        display('Rerunning uniquetol() & intersections.')
        % no stored copy or checksums were bad, must run uniquetol()s

        perpcs = newperpcs; % store checksums
        paracs = newparacs;

        [ v_vperpvals, v_vperpinds ] = uniquetol(v_perp,0.000001,'OutputAllIndices',true);
        [ v_vparavals, v_vparainds ] = uniquetol(v_para,0.000001,'OutputAllIndices',true);

        ss_uniques.v_vperpvals = v_vperpvals; ss_uniques.v_vperpinds = v_vperpinds;
        ss_uniques.v_vparavals = v_vparavals; ss_uniques.v_vparainds = v_vparainds;

        n_vperp = length(v_vperpvals);
        n_vpara = length(v_vparavals);

        % Make a grid for all possible (v_perp,v_para) tuples
        m_intersect = cell(n_vperp,n_vpara);
        m_interlen = zeros(n_vperp,n_vpara);
        for i=1:n_vperp
            parfor j=1:n_vpara

                % v_indices = intersect(v_vperpinds{i},v_vparainds{j});
                % using ismember() is faster, but still pretty slow
                m_intersect{i,j} = v_vperpinds{i}(ismember(v_vperpinds{i}, v_vparainds{j}));
                m_interlen(i,j) = length(m_intersect{i,j});

            end
        end

        % flatten
        m_intersect = reshape(m_intersect,1,[]);
        m_interlen = reshape(m_interlen,1,[]);

        % keep only points with matching cells
        m_intersect = m_intersect(m_interlen ~= 0);

        sm_intersect = m_intersect;
    end

    s_uniques = ss_uniques;
    m_intersect = sm_intersect;
    
end