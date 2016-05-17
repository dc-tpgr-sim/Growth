function [ maxw_vel ] = maxwellian( temp, shift_eV, PAcenter, PAwidth, in_velocities, in_angles)
%maxwellian(temperature, PA center, PA width, input eneriges, input angles) 
% Returns a discretely sampled, joint probability distribution function, 
% based on the input parameters, sampled at the provided energies and angles.
% Temp in K, shift in eV, PA in degrees, leave PAcenter empty [] to use a
% flat pitch-angle distribution.

    
    v_th = sqrt(3*temp*15156333.1);     % Convert input temp. to v_th = (3kT/m)^(1/2)
    eVconst = 3.913903e-6; % 2/(m_e*c^2) in eV^-1, i.e. conversion from eV to PSL
    v0 = 0.00989179273; % velocity base in PSL is equivalent to 25 eV

    shift = sqrt(shift_eV*eVconst)*299792458; % convert shift from eV to m/s

    % Maxwell-Boltzmann in velocity
    % maxw_vel = (temp/pi)^(3/2) * 4*pi * (in_velocities).^2 .* exp(-temp*(in_velocities-shift).^2);
    maxw_exp = exp(-(in_velocities-shift).^2/(2*v_th^2));
    maxw_vel = (2*pi)^(-3/2)*v_th^-3 .* maxw_exp;
    %sum(maxw_vel)
%    maxw_vel(maxw_vel<0) = 0;
    
    if ~isempty(PAcenter)
        % Gaussian in pitch angle
        maxw_PA = 1/(PAcenter*sqrt(2*pi)) * exp(-(in_angles-PAwidth).^2/(2*PAcenter^2));
        maxw_vel = maxw_vel.*maxw_PA;
    end
    
%    display([ 'shift: ' num2str(shift_eV) ' ds ' num2str(sum(maxw_vel)) ])
%    figure(fix(rand()*1000+50000))
%    plot(in_velocities, maxw_exp,'*')


end