% Redo for Run 4

%% Mirror data raw output manipulation, Run 2

Np = size(r_savX, 3);
B0 = 50e-6; % Magnetic field base is 50 uT
v0 = 0.00989179273; % likewise velocity base in PSL is equivalent to 25 eV
r0 = 0.337212985; % based on Larmour radius w/ above, length base is ~0.337 m
t0 = 7.14477319e-7; % based on B, Larmour period ~714 ns in s
dt = 0.01;
%% Gyro-interpolation

[ t_tz_time, t_Xf, t_Vf, t_mphi, t_circ ] = gyroterpolate(r_savX,r_savV,r_savB,0,dt,1:Np);

%% Build new results matrix

r_interp = zeros(3, 3, Np);

% new layout:
% 1,2 position and velocity at target-z
% 3,  number of timesteps before target (unitless, # of timesteps), 
    % fractional timestep to target (unitless, simulatin time [timesteps*dt]), 
    % total time (in ns, timesteps*dt*t0)

r_interp(:,1,:) = t_Xf.';
r_interp(:,2,:) = t_Vf.';

t_base_time = squeeze(r_res(2,3,:))-1;
t_time_ns = (t_base_time*dt+t_tz_time)*t0;

r_interp(:,3,:) = [ squeeze(r_res(2,3,:))-1 t_tz_time t_time_ns ].';

%% Hemi_fill to build a gyrotropic distribution.

t_dphi = 3*pi/256; % delta for co-latitude
t_domega = 0.001; % delta for solid angle in steradians
t_phis = 0+t_dphi:t_dphi:pi/2-t_dphi; % range of phis, discard first (pole) and last (plane)
    
% final data is stored flat in a 13xN matrix
% top x y z top vx vy vz bottom x y z bottom vx vy vz time (ns)
[ r_hemiterp, ~, r_hemirel ] = hemi_fill(r_interp, r_dist, t_dphi, t_domega);

%%

figure(8000)
t_en = sqrt(r_hemiterp(4,:).^2 + r_hemiterp(5,:).^2 + r_hemiterp(6,:).^2);
t_fr = sqrt(r_hemiterp(4,:).^2 + r_hemiterp(5,:).^2)./r_hemiterp(6,:);
scatter3(r_hemiterp(1,:),r_hemiterp(2,:),t_fr,[],t_en,'.')
%axis equal

%%

r_relfilt = find(~ismember(r_hemirel,t_zind));
r_hemifilt = r_hemiterp(:,r_relfilt);

figure(8001)
set(gcf,'position',[100 100 1000 600])
subplot(1,2,1)
plot(squeeze(sqrt(r_interp(1,2,:).^2+r_interp(2,2,:).^2)),squeeze(r_interp(3,2,:)),'.')
title('Full Set')
ylabel('$v_\parallel$','Interpreter','latex')
xlabel('$v_\perp$','Interpreter','latex')
subplot(1,2,2)
plot(squeeze(sqrt(r_interp(1,2,t_nzind).^2+r_interp(2,2,t_nzind).^2)),squeeze(r_interp(3,2,t_nzind)),'.')
title('Travel Time = 0 Particles Removed')
ylabel('$v_\parallel$','Interpreter','latex')
xlabel('$v_\perp$','Interpreter','latex')

scatter(t_angles,t_energies)
title('Energy vs. Pitch Angle for Travel Time = 0 Particles')
ylabel('Energy [eV]'); set(gca,'ytick',unique(t_energies));
xlabel('Pitch Angle [deg]'); 


%% Save

save J:\Data' Core'\particles\sim\run4-hemi.mat t_domega t_dphi t_phis ...
    B0 r0 t0 v0 target_length target_z N_part dt mirror_ratio length_factor ...
    r_dist r_res r_interp r_hemiterp 

%% Convert to units and reverse all velocities

Np = size(r_hemiterp,2);
eVconst = 3.913903e-6;

t_X_t = r_hemiterp(1:3,:)*r0; % Distances in m

t_vx_t = -r_hemiterp(4,:);
t_vy_t = -r_hemiterp(5,:);
t_vz_t = -r_hemiterp(6,:);
t_vmag_t = sqrt(t_vx_t.^2 + t_vy_t.^2 + t_vz_t.^2); % v, unitless
t_vpar_t = t_vz_t; % vpar, unitless
t_vper_t = sqrt(t_vx_t.^2 + t_vy_t.^2); % vper, unitless

t_alpha_t = atan2(t_vper_t,-t_vpar_t);%*180/pi;
t_theta_t = 2*pi-atan2(t_vy_t,t_vx_t); % math to convert from atan2 output range to standard
t_theta_t(t_theta_t>=2*pi) = t_theta_t(t_theta_t>=2*pi)-2*pi; % 0to2pi clockwise from +x angle
%t_theta_t = t_theta_t * 180/pi;

t_X_b = r_hemiterp(7:9,:)*r0; % Distances in m

t_En_t = (t_vmag_t*v0).^2/eVconst; % Energy, eV
t_vmag_t_mps = t_vmag_t*v0*299792458; % convert to m/s
t_vpar_t_mps = t_vpar_t*v0*299792458;
t_vper_t_mps = t_vper_t*v0*299792458;

t_vx_b = -r_hemiterp(10,:);
t_vy_b = -r_hemiterp(11,:);
t_vz_b = -r_hemiterp(12,:);
t_vmag_b = sqrt(t_vx_b.^2 + t_vy_b.^2 + t_vz_b.^2); % v, unitless
t_vpar_b = t_vz_b; % vpar, unitless
t_vper_b = sqrt(t_vx_b.^2 + t_vy_b.^2); % vper, unitless

t_alpha_b = atan2(t_vper_b,-t_vpar_b);%*180/pi;
t_theta_b = 2*pi-atan2(t_vy_b,t_vx_b); % math to convert from atan2 output range to standard
t_theta_b(t_theta_b>=2*pi) = t_theta_b(t_theta_b>=2*pi)-2*pi; % 0to2pi clockwise from +x angle
%t_theta_b = t_theta_b * 180/pi;

t_En_b = (t_vmag_t*v0).^2/eVconst; % Energy, eV
t_vmag_b_mps = t_vmag_b*v0*299792458; % convert to m/s
t_vpar_b_mps = t_vpar_b*v0*299792458;
t_vper_b_mps = t_vper_b*v0*299792458;

% 17xN matrix
% (x,y,z, vmag, vpar, vper, pa, azi)_top
% (x,y,z, vmag, vpar, vper, pa, azi)_bottom
% time
r_mirror_XVPA = [ ... 
    t_X_t(1,:) ; t_X_t(2,:) ; t_X_t(3,:) ; ...
    t_vmag_t_mps ; t_vpar_t_mps ; t_vper_t_mps ; ...
    t_alpha_t ; t_theta_t ; ...
    t_X_b(1,:) ; t_X_b(2,:) ; t_X_b(3,:) ; ...
    t_vmag_b_mps ; t_vpar_b_mps ; t_vper_b_mps ; ...
    t_alpha_b ; t_theta_b ; ...
    r_hemiterp(13,:) ...
    ];

% Discard X (assume homogeneity), readd Energies, 13xN matrix
% 13xN matrix
% (En, vmag, vpar, vper, pa, azi)_top
% (En, vmag, vpar, vper, pa, azi)_bottom
% time
r_mirror_EVPA = [ ... 
    t_En_t ; ...
    t_vmag_t_mps ; t_vpar_t_mps ; t_vper_t_mps ; ...
    t_alpha_t ; t_theta_t ; ...
    t_En_b ; ...
    t_vmag_b_mps ; t_vpar_b_mps ; t_vper_b_mps ; ...
    t_alpha_b ; t_theta_b ; ...
    r_hemiterp(13,:) ...
    ];

smap_EVPA.top.En = 1;
smap_EVPA.top.v.mag = 2;
smap_EVPA.top.v.para = 3;
smap_EVPA.top.v.perp = 4;
smap_EVPA.top.alpha = 5;
smap_EVPA.top.theta = 6;
smap_EVPA.bot.En = 7;
smap_EVPA.bot.v.mag = 8;
smap_EVPA.bot.v.para = 9;
smap_EVPA.bot.v.perp = 10;
smap_EVPA.bot.alpha = 11;
smap_EVPA.bot.theta = 12;
smap_EVPA.time = 13;

