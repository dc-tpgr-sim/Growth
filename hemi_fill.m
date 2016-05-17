function [ n, n_azi, rel ] = hemi_fill(xvt,r_dist,t_dphi,t_domega,varargin)
% hemi_fill: function to populate a constant-solid-angle hemisphere, given
% a single stripe of co-latitude positions and velocities, the
% corresponding co-latitudes, and the solid angle value in steradians.
%
% xv should be a 3x2xN vector, where N=(2pi/dphi)-2, 
% positions are in (:,1,:), and velocities in (:,2,:)

    opt = struct('cell',false,'phistop',2*pi);
    opt = optParse(opt,varargin{:});

    % range of thetas, discard first (pole) and last (equator)
    t_phis = 0+t_dphi:t_dphi:pi/2-t_dphi;
    n_phis = length(t_phis);
    t_alpha = sqrt(r_dist(4,:).^2 + r_dist(5,:).^2)./r_dist(6,:);
    n_En = length(find(t_alpha == 0));
    
    nc_Xt = cell(n_phis,n_En,1);
    nc_Vt = cell(n_phis,n_En,1);
    nc_Xb = cell(n_phis,n_En,1);
    nc_Vb = cell(n_phis,n_En,1);
    nc_t = cell(n_phis,n_En,1);
    n_azi = zeros(n_phis,1);
    nc_rel = cell(n_phis,n_En,1);
    for i=1:n_phis
        t_phi = t_phis(i);

        n_azi(i) = round(2*pi*sin(t_phi)*t_dphi/t_domega);
        l_thetas = 0:2*pi/n_azi(i):opt.phistop;
        n_az = length(l_thetas);
        [ n_az n_azi(i) ];
        if n_az ~= n_azi(i)
        %    display('fuuu');
            n_azi(i) = n_az;
        end
        for j=1:n_En
            part = (i-1)*n_En + j;
            display(['fnh ' num2str(part) ' lsjdf ' num2str(size(xvt))])
            t_x = squeeze(xvt(:,1,part));
            t_v = squeeze(xvt(:,2,part));
            b_x = r_dist(1:3,part);
            b_v = r_dist(4:6,part);
            t_t = squeeze(xvt(:,3,part));        
    
            n_xt = zeros(3,n_az);
            n_vt = zeros(3,n_az);
            n_xb = zeros(3,n_az);
            n_vb = zeros(3,n_az);
            for k=1:n_az
            
                t_theta = l_thetas(k);
                t_rot = [ cos(t_theta) sin(t_theta) 0 ; -sin(t_theta) cos(t_theta) 0 ; 0 0 1 ];
                
                n_xt(:,k) = t_rot*t_x;
                n_vt(:,k) = t_rot*t_v;
                n_xb(:,k) = t_rot*b_x;
                n_vb(:,k) = t_rot*b_v;
            end
            
%            n_xt
%            n_vt
%            n_xb
%            n_vb
%            t_t
            
            nc_Xt{i,j} = n_xt;
            nc_Vt{i,j} = n_vt;
            nc_Xb{i,j} = n_xb;
            nc_Vb{i,j} = n_vb;
            nc_t{i,j} = t_t(3);
            nc_rel{i,j} = part;
        end
                
    end
    
    n_vec = sum(n_azi)*n_En;
    if opt.cell
        n = cell(n_phis,5);
        n(:,1) = nc_Xt;
        n(:,2) = nc_Vt;
        n(:,3) = nc_Xb;
        n(:,4) = nc_Vb;
        n(:,5) = nc_t;
    else
        n = zeros(13,n_vec);
        rel = zeros(1,n_vec);
        i_n = 0;
        for i=1:n_phis
            for j=1:n_En
                for k=1:n_azi(i)
                    i_n = i_n + 1;
%                    size([ nc_Xt{i,j}(:,k) ; nc_Vt{i,j}(:,k) ; nc_Xb{i,j}(:,k) ; nc_Vb{i,j}(:,k) ; nc_t{i,j} ])
                	n(:,i_n) = [ nc_Xt{i,j}(:,k) ; nc_Vt{i,j}(:,k) ; nc_Xb{i,j}(:,k) ; nc_Vb{i,j}(:,k) ; nc_t{i,j} ];
                    rel(i_n) = nc_rel{i,j};
                end
            end
        end
    end
    
end

function optstr = optParse(options, varargin)

    %# read the acceptable names
    optionNames = fieldnames(options);

    %# count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('hemi_fill needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
       inpName = lower(pair{1}); %# make case insensitive

       if any(strcmp(inpName,optionNames))
          %# overwrite options. If you want you can test for the right class here
          %# Also, if you find out that there is an option you keep getting wrong,
          %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    
    optstr = options;
end