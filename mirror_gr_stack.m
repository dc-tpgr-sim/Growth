% Full distribution to growth rate code stack

%% Define the distribution

% Fiddle with topside distribution
fiddle = false;
if fiddle
    % Find unique energies for fiddling
    [ t_en, t_en_ind ] = uniquetol(r_mirror_EVPA(smap_EVPA.top.En,:),0.001,'outputallindices',true);
    
    % Remove half of the energies
    t_en_find = vertcat( t_en_ind{1:2:end} );
    GR_dist = r_mirror_EVPA(:,t_en_find);
else
    GR_dist = r_mirror_EVPA;
end

GR_smap = smap_EVPA;

shortest_travel_time = min(GR_dist(GR_smap.time,:));
longest_travel_time = max(GR_dist(GR_smap.time,:));

density_const = 0.000314207783; % m_e*epsilon_0/e^2

% ionospheric background parameters
f_pe = 400000;
omega_pe = f_pe * 2*pi;
n_e = omega_pe^2*density_const;
iono_temperature = 2000; %in Kelvin

s_dyn_iono = struct('n', n_e, 'temp', iono_temperature, 'shift', 0, ...
    'PAcenter', [], 'PAwidth', []);

% secondary background parameters
maxw_temperature = 200000; % in Kelvin
t_temp_eV = maxw_temperature/11604.505;
if t_temp_eV < 137.1 % linearly interpolate n(T) from Table 1b in Lotko & Maggs 1981
    t_n_lotkomaggs = (0.44 - 0.83)/(137.1 - 62.4)*(t_temp_eV-62.4) + 0.83;
else
    t_n_lotkomaggs = (0.42 - 0.44)/(220.1 - 137.1)*(t_temp_eV-137.1) + 0.44;
end
maxw_particles = t_n_lotkomaggs * 1000000; % cm^-3 -> m^-3
maxw_shift = 0;
maxw_PAcenter = [];
maxw_PAwidth = [];

s_dyn_bg = struct('n', maxw_particles, 'temp', maxw_temperature, 'shift', maxw_shift, ...
    'PAcenter', maxw_PAcenter, 'PAwidth', maxw_PAwidth);

% beam parameter sets
s_dyn_beams = {
    struct('n', 0, 'dwell_time', 5), ... % run 4 longest travel time is <14s
    struct('n', maxw_particles/50, 'dwell_time', 5, 'temp', maxw_temperature/5, 'shift', 400, 'PAcenter', [], 'PAwidth', []), ...
...%    struct('n', maxw_particles/5, 'temp', maxw_temp/16, 'shift', 300, 'PAcenter', [], 'PAwidth', []), ...
    struct('n', 0, 'dwell_time', 5)
    };

% Builder function eats launch time, sample time, and the three dist structures.
s_dyn_dist = build_dyn_struct(0.001, 0.010, s_dyn_iono, s_dyn_bg, s_dyn_beams);

%% Sanity-check plots

java_numFmt = java.text.DecimalFormat;

% sort by energy for plotting, but we feed dynamic_distribution velocities and PAs
[ s_En, si_En ] = sort(GR_dist(GR_smap.bot.En,:));

h = figure(7777);
clf(h)

set(h,'position',[ 10 500 300*n_beams 400]);
suptitle('Electron Distribution Functions')

n_beams = length(s_dyn_beams);
for i=1:n_beams
        
    [ t_dist, s_dist ] = dynamic_distribution(s_dyn_dist.times(i), GR_dist, GR_smap, s_dyn_dist);

    t_width = 0.90/n_beams;
    subplot('position',[ 0.05+t_width*(i-1) 0.17 t_width 0.70 ])

    plot(s_En,t_dist(si_En),'k',s_En,s_dist{1}(si_En),'g.',s_En,s_dist{2}(si_En),'b*',s_En,s_dist{3}(si_En),'rx');
    set(gca,'fontsize',12)
    
    xlabel('Energy [eV]')
    foo = get(gca,'xticklabel'); foo{end}=''; set(gca,'xticklabel',foo);
    if i==1
        ylabel('$f(|v|)$','interpreter','latex');
        t_xlim = xlim; t_ylim = ylim;
    else
        set(gca,'yticklabel',[])
        xlim(t_xlim); ylim(t_ylim);
    end
    
    if s_dyn_dist.beams{i}.n == 0
        legend('Combined', [ char(java_numFmt.format(s_dyn_dist.iono.temp)) ' K ionospheric BG' ], ...
            [ char(java_numFmt.format(s_dyn_dist.bg.temp)) ' K secondary BG' ])
    else
        legend('Combined', [ char(java_numFmt.format(s_dyn_dist.iono.temp)) ' K ionospheric BG' ], ...
            [ char(java_numFmt.format(s_dyn_dist.bg.temp)) ' K secondary BG' ], ...
            [ char(java_numFmt.format(s_dyn_dist.beams{i}.temp)) ' K, ' ...
                char(java_numFmt.format(s_dyn_dist.beams{i}.shift)) ' eV-shifted beam'])
    end
end

%print('-dpng',[file_outdir '\topdist.png'])
%close(h)


%% azi_sum all top timesteps
launch_time = s_dyn_dist.launch_dt;
n_beams = length(s_dyn_dist.beams);

tic
dt_c_EVPN = cell(n_launchsteps,2);

v_launchsteps = 0:launch_time:s_dyn_dist.times(end)-launch_time;
n_launchsteps = length(v_launchsteps);

for i=1:n_launchsteps
    t_time = v_launchsteps(i);
    
    [ m_EVPN, smap_EVPN ] = time_azi_sum_chain (t_time, GR_dist, GR_smap, s_dyn_dist);
    
    t_strike = m_EVPN(smap_EVPN.time, :) + t_time;
    
    dt_c_EVPN{i,1} = m_EVPN;
    dt_c_EVPN{i,2} = t_strike;
    
end

dt_m_EVPN = [dt_c_EVPN{:,1}];
dt_v_EVPNt = [dt_c_EVPN{:,2}];
toc

%% Set up for perp_sum
% Now we go through and filter, for bottom time t-deltat to t

sample_time = s_dyn_dist.sample_dt;
v_timesteps = shortest_travel_time:sample_time:s_dyn_dist.times(end)-launch_time;
n_timesteps = length(v_timesteps);

% Find the field-aligned velocities, for use as the center points in the
% reduced distribution function.
display('Uniquetol...')
v_paravels = uniquetol(dt_m_EVPN(smap_EVPN.bot.v.para,dt_m_EVPN(smap_EVPN.bot.v.perp,:)==0));
display('...done.')

% Reverse to smallest magnitude first
v_paravels = sortmag(v_paravels);

% extend these to zero
n_para = length(v_paravels);
d_vpar = median(diff(v_paravels));
d_extrap = v_paravels(1):-d_vpar:0;

v_paravelx = [ flip(d_extrap(2:end)) v_paravels ];
n_paravelx = length(v_paravelx);

%% perp_sum all top timesteps

display('Running perp_sum()s...')
tic
dt_m_rdf = zeros(n_timesteps,n_paravelx);
dt_v_nrdf = zeros(n_timesteps,1);
for i=1:n_timesteps
    
    t_time = v_timesteps(i);
    % search for particles that have been 'detected' in this timeslice
    v_timeinds = find(dt_v_EVPNt <= t_time & dt_v_EVPNt > t_time-sample_time);    
    m_particles = dt_m_EVPN(:,v_timeinds);
    dt_v_nrdf(i) = size(m_particles,2);
    m_particles(smap_EVPN.dist,:) = m_particles(smap_EVPN.dist,:) * launch_time/sample_time;
    
    if length(m_particles) < 1
        continue
    end
    
%    display(['Time ' num2str(t_time) ' found ' num2str(length(m_particles)) ' particles.'])
    
    % reduce to parallel
    dt_m_rdf(i,:) = perp_sum(m_particles, smap_EVPN, v_paravelx);

    % growth rate
%    v_gRate = para_gRate(m_rdf);
    
end
toc

%%

for i=490:800 %1:n_timesteps
    
    t_time = v_timesteps(i);
    % search for particles that have been 'detected' in this timeslice
    [ t_time t_time-sample_time ]
    v_timeinds = find(dt_v_EVPNt <= t_time & dt_v_EVPNt > t_time-sample_time);
    length(v_timeinds)
end

%% Full gamma vs k & time plot

t_temp = s_dyn_dist.beams{2}.temp;
t_shift = s_dyn_dist.beams{2}.shift;
t_bg = 2000; % background ionospheric cold electron temperature [K]
v_bg = sqrt(3*t_bg*15156333.1); 

f_pe = 400000; % 500 kHz plasma freq.
omega_pe = f_pe * 2*pi;
t_test_temp = t_shift + t_temp/11604; % approximate beam speed
t_test_vel = eV2mps(t_test_temp);
[ ~, i_test ] = min(abs(-v_paravelx-t_test_vel))
v_omega_test = (1.00001:0.00001:1.01)*omega_pe;
v_k_test = sqrt(2/3*(v_omega_test-omega_pe)*omega_pe/v_bg^2);
%v_k_test = logspace(-6,20,1000);
v_k_test = 0.1:0.001:0.5;
v_omega_test = v_k_test.^2*3/2*v_bg^2/omega_pe + omega_pe;
n_test = length(v_k_test)

m_gamma = zeros(n_timesteps,n_test); % timestep,kind,val/omegaind
m_vtest = zeros(n_timesteps,n_test,2);
m_kmag2 = zeros(n_timesteps,n_test,1);
m_df1 = zeros(n_timesteps,n_test,1);
m_omega_test = zeros(n_timesteps,n_test,1);
m_n_e = zeros(n_timesteps,n_test,1);

parfor i=1:n_timesteps
    for j=1:n_test
        t_kpara = v_k_test(j);
%        t_omega_test = v_omega_test(j);

        [ m_gamma(i,j), m_vtest(i,j,:), m_kmag2(i,j), m_df1(i,j), m_omega_test(i,j), m_n_e(i,j) ] = growth_rate(dt_m_rdf(i,:), -v_paravelx, [ t_kpara 0 ], omega_pe, v_bg, t_test_vel);

    end
end

%%

subplot(2,2,1)
plot(v_k_test,m_kmag2(750,:)); title('kmag2')
subplot(2,2,2)
plot(v_k_test,m_df1(750,:)); title('df1')
subplot(2,2,3)
plot(v_k_test,m_omega_test(750,:)); title('omega test')
subplot(2,2,4)
plot(v_k_test,m_n_e(750,:)); title('n_e')

%% Launch timestep n summation, for 'this is where the beam was' plot.

dt_v_fbg = zeros(n_launchsteps,1);
dt_v_fbeam = zeros(n_launchsteps,1);
parfor i=1:n_launchsteps
    t_time = v_launchsteps(i);

    [ ~, s_dist ] = dynamic_distribution(t_time, GR_dist, GR_smap, s_dyn_dist);

    dt_v_fbg(i) = sum(s_dist{1}) + sum(s_dist{2});
    dt_v_fbeam(i) = sum(s_dist{3});
end


%% r/b gamma vs k,time plot

t_azi = 0;
t_el = 0;

v_upsteps = find(v_timesteps > 7.2 & v_timesteps < 8);
v_dnsteps = find(v_timesteps >= 6 & v_timesteps < 12);

%v_upsteps = find(v_timesteps);
%v_dnsteps = find(v_timesteps);

h = figure(7805);
clf
set(h, 'position', [100 50 1200 900])

p_plbase = 0.08;
p_plleft = 0.08;
p_plwidt = 0.40;
p_lpheig = 0.10;
p_vspace = 0.03;
p_spheig = 0.32;
p_hspace = 0.04;

hT = suptitle([ 'Growth Rates, $\Delta t_S =' num2str(launch_time) '$ s, $\Delta t_D =' num2str(sample_time) '$ s' ]);
set(hT,'interpreter','latex');

% -- n vs t --

nax = subplot('Position',[ ...
    p_plleft ...
    p_plbase+2*p_spheig+p_vspace ...
    2*p_plwidt+p_hspace ...
    p_lpheig ...
]); ax = [ ax nax ];
plot(v_launchsteps,dt_v_fbeam./dt_v_fbg)
xlabel('Time [s]'); set(gca, 'fontsize', 12); grid on; ylim([-0.25 0.75]); set(gca,'ytick',[0 0.2 0.4 0.6])
set(gca,'XAxisLocation','top'); ylabel('$n_{beam}/n_{bg}$','interpreter','latex')
set(gca,'xtick',1:15)

% -- gamma vs t --

ax = [];
nax = subplot('Position',[ ...
    p_plleft ...
    p_plbase+p_spheig ...
    p_plwidt ...
    p_spheig ...
]); ax = [ ax nax ];
surf(v_timesteps(v_upsteps),v_k_test,m_gamma(v_upsteps,:).','edgecolor','none'); colormap(rwbmap); box on;  set(gca, 'layer', 'top')
%caxis([-t_crange t_crange])
view(0,0); %set(gca,'zscale', 'log');  % ylim([min(v_k_test) 0.5])
%zlim([-10e3 10e3])
set(gca, 'fontsize', 12, 'xticklabel', []); ylabel('k'); zlabel('\gamma');
t_tick = get(gca,'ztick'); t_tick = t_tick(2:end); set(gca, 'ztick',t_tick);

nax = subplot('Position',[ ...
    p_plleft+p_hspace+p_plwidt ...
    p_plbase+p_spheig ...
    p_plwidt ...
    p_spheig ...
]); ax = [ ax nax ];
surf(v_timesteps(v_dnsteps),v_k_test,m_gamma(v_dnsteps,:).','edgecolor','none'); colormap(rwbmap); box on; set(gca, 'layer', 'top')
t_crange = max([caxis(ax(1)) caxis(ax(2))]);
caxis(ax(1),[-t_crange t_crange]); caxis(ax(2),[-t_crange t_crange])
view(0,0); set(gca, 'ydir', 'reverse'); %set(gca, 'zscale', 'log');% ylim([min(v_k_test) 0.5])
zlim([-10e3 10e3])
set(gca, 'fontsize', 12, 'xticklabel', []); %ylabel('k'); zlabel('\gamma'); 
t_tick = get(gca,'ztick'); t_tick = t_tick(2:end); set(gca, 'ztick',t_tick);

% -- k vs t --

nax = subplot('Position',[ ...
    p_plleft ...
    p_plbase ...
    p_plwidt ...
    p_spheig ...
]); ax = [ ax nax ];
surf(v_timesteps(v_upsteps),v_k_test,m_gamma(v_upsteps,:).','edgecolor','none'); colormap(rwbmap); box on; set(gca, 'layer', 'top')
%caxis([-t_crange t_crange])
view(0,90); ylim([min(v_k_test) 0.5])
%ylim([0 3e-3])
set(gca, 'fontsize', 12,'xtick',get(ax(1),'xtick'))
ylabel('k');
t_tick = get(gca,'yticklabel'); t_tick{end}=''; set(gca, 'yticklabel',t_tick);
xlabel('Time [s]');

nax = subplot('Position',[ ...
    p_plleft+p_hspace+p_plwidt ...
    p_plbase ...
    p_plwidt ...
    p_spheig ...
]); ax = [ ax nax ];
surf(v_timesteps(v_dnsteps),v_k_test,m_gamma(v_dnsteps,:).','edgecolor','none'); colormap(rwbmap); box on; set(gca, 'layer', 'top')
t_crange = max([caxis(ax(3)) caxis(ax(4))]);
caxis(ax(3),[-t_crange t_crange])
caxis(ax(4),[-t_crange t_crange])
view(0,-90); set(gca, 'ydir', 'reverse'); ylim([min(v_k_test) 0.5])
%ylim([0 3e-3])
set(gca, 'fontsize', 12);  %ylabel('k')
t_tick = get(gca,'yticklabel'); t_tick{end}=''; set(gca, 'yticklabel',t_tick);
xlabel('Time [s]');

foo = annotation('line',[0.08 0.473],[0.72 0.749]);
set(foo,'color','red')
uistack(foo,'bottom')
foo = annotation('line',[0.48 0.528],[0.72 0.749]);
set(foo,'color','red')
uistack(foo,'bottom')

foo = annotation('line',[0.521 0.752],[0.72 0.749]);
set(foo,'color','red');
uistack(foo,'bottom')
foo = annotation('line',[0.92 0.808],[0.72 0.749]);
set(foo,'color','red');
uistack(foo,'bottom')

foo = annotation('textbox',[0.45 0.325 0.1 0.1], ...
    'string', '...', 'horizontalalignment', 'center', 'linestyle', 'none','fontsize', 16, 'fontweight', 'bold');

% print looks awful, use manual export
%print('-opengl','-dpng', [file_outdir '\gr.png'])


%% r/b gamma vs k,time plot, single zoom

t_azi = 0;
t_el = 0;

v_steps = find(v_timesteps > 7.2 & v_timesteps < 8);

h = figure(7805);
clf
set(h, 'position', [100 50 1200 900])

p_plbase = 0.08;
p_plleft = 0.08;
p_plwidt = 0.84;
p_lpheig = 0.10;
p_vspace = 0.03;
p_spheig = 0.32;

hT = suptitle([ 'Growth Rates, $\Delta t_S =' num2str(launch_time) '$ s, $\Delta t_D =' num2str(sample_time) '$ s, 100 ms Beam' ]);
set(hT,'interpreter','latex');

% -- n vs t --

nax = subplot('Position',[ ...
    p_plleft ...
    p_plbase+2*p_spheig+p_vspace ...
    p_plwidt ...
    p_lpheig ...
]); ax = [ ax nax ];
plot(v_launchsteps,dt_v_fbeam./dt_v_fbg)
xlabel('Time [s]'); set(gca, 'fontsize', 12); grid on; ylim([-0.25 0.75]); set(gca,'ytick',[0 0.2 0.4 0.6])
set(gca,'XAxisLocation','top'); ylabel('$n_{beam}/n_{bg}$','interpreter','latex')
set(gca,'xtick',[ 1:7 8:10]); xlim(v_launchsteps([1 end]));

% -- gamma vs t --

ax = [];
nax = subplot('Position',[ ...
    p_plleft ...
    p_plbase+p_spheig ...
    p_plwidt ...
    p_spheig ...
]); ax = [ ax nax ];
surf(v_timesteps(v_upsteps),v_k_test,m_gamma(v_upsteps,:).','edgecolor','none'); colormap(rwbmap); box on;  set(gca, 'layer', 'top')
t_crange = max([abs(caxis(ax(1)))]);
caxis([-t_crange t_crange])
view(0,0); %set(gca,'zscale', 'log');  % ylim([min(v_k_test) 0.5])
%zlim([-10e3 10e3])
set(gca, 'fontsize', 12, 'xticklabel', []); ylabel('k'); zlabel('\gamma');
t_tick = get(gca,'ztick'); t_tick = t_tick(2:end); set(gca, 'ztick',t_tick);

% -- k vs t --

nax = subplot('Position',[ ...
    p_plleft ...
    p_plbase ...
    p_plwidt ...
    p_spheig ...
]); ax = [ ax nax ];
surf(v_timesteps(v_upsteps),v_k_test,m_gamma(v_upsteps,:).','edgecolor','none'); colormap(rwbmap); box on; set(gca, 'layer', 'top')
%t_crange = max([abs(caxis(ax(2)))]);
caxis([-t_crange t_crange])
view(0,90); ylim([min(v_k_test) 0.5])
%ylim([0 3e-3])
set(gca, 'fontsize', 12,'xtick',get(ax(1),'xtick'))
ylabel('k');
t_tick = get(gca,'yticklabel'); t_tick{end}=''; set(gca, 'yticklabel',t_tick);
xlabel('Time [s]');

foo = annotation('line',[0.08 0.6805],[0.72 0.749]);
set(foo,'color','red')
uistack(foo,'bottom')
foo = annotation('line',[0.92 0.745],[0.72 0.749]);
set(foo,'color','red')
uistack(foo,'bottom')

foo = annotation('textbox',[0.45 0.325 0.1 0.1], ...
    'string', '...', 'horizontalalignment', 'center', 'linestyle', 'none','fontsize', 16, 'fontweight', 'bold');

% print looks awful, use manual export
%print('-opengl','-dpng', [file_outdir '\gr.png'])
