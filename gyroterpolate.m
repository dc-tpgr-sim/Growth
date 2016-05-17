function [ t_tz_time, t_Xf, t_Vf, t_mphi, t_circ ] = gyroterpolate(t_X, t_V, t_B, target_z, dt, t_d)
% Takes 3xTxN input vectors, where T = some number of saved timesteps, and
% N = some number of particles which have crossed the target z-level.
% Fits a gyroorbit to each particle's track, then interpolates the actual
% strike XVT.

    Ns = size(t_X, 2);
    Np = size(t_X, 3);
    t_tz_time = zeros(Np,1);
    t_Xf = zeros(Np,3);
    t_Vf = zeros(Np,3);
    t_mphi = zeros(Np,1);
    t_circ = zeros(Np,3);
    
    options = optimoptions('fminimax');
    options.Display = 'none';
    
    parfor part=t_d%1:Np
        
        t_pX = squeeze(t_X(:,:,part));
        t_pV = squeeze(t_V(:,:,part));
        
        tz_center = HyperSVD(squeeze(t_pX(1:2,:)).');
        
        if tz_center(3) ~= 0 % gyropath
            
            t_circ(part,:) = tz_center;

            % We know the equations of motion that the particle should be
            % following; the only thing we don't know is the phase.
            omega = sqrt(sum(squeeze(t_B(:,:,part)).^2,1)); % Bmag
            vperp = sqrt( squeeze(t_pV(1,:)).^2 + squeeze(t_pV(2,:)).^2 );
            vpar = t_pV(3,:);
            ttime = (-(Ns-2):1)*dt;

            t_Xc = [ t_pX(1,:) ; t_pX(2,:) ; t_pX(3,:) ];
            t_x0 = t_pX(1,end-1);
            t_y0 = t_pX(2,end-1);
            t_z0 = t_pX(3,end-1);
            deltx = t_Xc(1,:);
            delty = t_Xc(2,:);
            deltz = t_Xc(3,:)-t_z0;
            Cx = -vperp./omega;
            Cy =  vperp./omega;
            Cz = vpar.*ttime;
            tau = omega.*ttime;

            ferrorphi = @(phi) ferrorphifunc(phi,deltx,delty,deltz,tau,Cx,Cy,Cz); % geometric error

            [ t_mphi(part), ~ ] = fminimax(ferrorphi,pi,[],[],[],[],0,2*pi,[],options);

            tz_x0 = t_pX(1,end-1); tz_vx0 = t_pV(1,end-1); % X, V before target crossing
            tz_y0 = t_pX(2,end-1); tz_vy0 = t_pV(2,end-1);
            tz_z0 = t_pX(3,end-1); tz_vz0 = t_pV(3,end-1);
            t_tz_time(part) = (target_z-tz_z0)/vpar(end-1); % travel time to target crossing

            % now just use gryo equations to get final interp. results
            tz_time = t_tz_time(part);
            tz_vperp = sqrt(tz_vx0^2 + tz_vy0^2);
            tz_omega = omega(end-1);
            tz_tau = tz_omega*(ttime(end-1)+tz_time);
            tz_phi = t_mphi(part);

            tz_xf = -tz_vperp/tz_omega*cos(tz_tau + tz_phi) + tz_center(1);
            tz_yf =  tz_vperp/tz_omega*sin(tz_tau + tz_phi) + tz_center(2);
            tz_zf =  tz_vz0*tz_time + tz_z0;

            tz_vxf = tz_vperp*sin(tz_tau + tz_phi);
            tz_vyf = tz_vperp*cos(tz_tau + tz_phi);
            tz_vzf = tz_vz0;

            [ tz_x0 tz_y0 tz_z0 ];
            t_Xf(part,:) = [ tz_xf tz_yf tz_zf ];
            [ tz_vx0 tz_vy0 tz_vz0 ];
            t_Vf(part,:) = [ tz_vxf tz_vyf tz_vzf ];
            
        else % straight line
            
            t_tz_time(part) = (target_z-t_pX(3,end-1))/t_pV(3,end-1); % travel time to target crossing
            
            % x, y, and velocities don't change, just set z = target_z
            t_Xf(part,:) = [ t_pX(1,end-1) t_pX(2,end-1) target_z ].';
            t_Vf(part,:) = t_pV(:,end-1).';
            t_mphi(part) = NaN;
            t_circ(part,:) = [ t_pX(1,end-1) t_pX(2,end-1) 0 ];
            
        end
            

        
    end
        
    %fdx = fdx + t_circ(1);
    %fdy = fdy + t_circ(2);
    
    %    phi = fminbnd(ferrorphi, 0, 2*pi);

end

function err = ferrorphifunc(phi,deltx,delty,deltz,tau,Cx,Cy,Cz)

    err = sqrt( ...
        (deltx - Cx.*cos(tau + phi)).^2 + ...
        (delty - Cy.*sin(tau + phi)).^2 + ...
        (deltz - Cz).^2);

end

function vp = vel_upd(t, v, B, phi)
    v_perp = sqrt(v(1).^2 + v(2).^2);
    omega_g = 2*pi*B;

    vp(1) = v_perp.*sin(omega_g.*t + phi);
    vp(2) = v_perp.*sin(omega_g.*t + phi);
    vp(3) = v(3);
end

function xp = pos_upd(t, x, v, B, phi)
    v_perp = sqrt(v(1).^2 + v(2).^2);
    omega_g = 2*pi*B;
    
    xp(1) = x(1) + -v_perp/omega_g.*cos(omega_g.*t + phi);
    xp(2) = x(2) +  v_perp/omega_g.*sin(omega_g.*t + phi);
    xp(3) = x(3) + v(3).*t;
end