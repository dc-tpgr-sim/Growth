function [ dist, sdist, n_beam ] = dynamic_distribution(time, in_dist, in_map, dist_def)
% Returns a distribution at a given time, for a provided 2xN list of 
% particles/distribution function element centers, with velocities in (1,:) 
% and pitch angles in (2,:).  Returns an N-element list which is values 
% of f(vmag,pa) for each particle.

    in_vmag = in_dist(in_map.bot.v.mag,:);
    in_PA = in_dist(in_map.bot.alpha,:);
    
    % ionosphere parameters
    iono_def = dist_def.iono;
    
    % Create ionospheric distribution
    iono_dist = maxwellian(iono_def.temp, iono_def.shift, ...
        iono_def.PAcenter, iono_def.PAwidth, in_vmag, in_PA);
    iono_part = iono_dist*iono_def.n;
    
    % background parameters
    bg_def = dist_def.bg;

    % Create background distribution
    bg_dist = maxwellian(bg_def.temp, bg_def.shift, ...
        bg_def.PAcenter, bg_def.PAwidth, in_vmag, in_PA);
    bg_part = bg_dist*bg_def.n;
    
%    display(sum(num2str(bg_dist)))
%    display(['bg n ' num2str(bg_def.n) ' ds ' num2str(sum(bg_part))])
    
    % We'll be using segment time /. dwell_time
    i_beam = find(dist_def.times <= time, 1, 'last');
    beam_def = dist_def.beams{i_beam};
    n_beam = beam_def.n;
        
    if beam_def.n == 0 % BG-only case

        dist = bg_part+iono_part;
        sdist = { iono_part, bg_part, zeros(size(iono_part)) };
        
    else % BG + beam
        
        beam_dist = maxwellian(beam_def.temp, beam_def.shift, ...
            beam_def.PAcenter, beam_def.PAwidth, in_vmag, in_PA);
        beam_part = beam_dist*beam_def.n;
%        display(sum(num2str(beam_dist)))
%        display(['beam n ' num2str(beam_def.n) ' ds ' num2str(sum(beam_part))])
        
        dist = iono_part+bg_part+beam_part;
        sdist = { iono_part, bg_part, beam_part };
        
    end

end