function  genMultipleAssociationResults_v2()
clc;close all; clear;
%format long e;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.P = 1;
params.simulation_area_side = [-100 100];    % simulation area side to side
params.space_realizations = 10;
params.time_slots = 1;
params.M = 1;
params.M_agg = 4;
sim = 'InterferersDistribution';
%sim ='CellsDensityADR';
%sim ='CellsDensityASE';
%sim ='CellsDensityASE_Aggregate';
%sim ='UsersDensityADR';
%sim ='UsersDensityASE';
%sim ='UsersDensityASE_Aggregate';
%sim ='CompareToRef';
%sim ='ImpactOfFading';
%sim ='AggregateRate';
%sim = 'OvelappedMultiCells';
%sim = 'RangeOfDistanceToNthNeighbour';
%sim = 'CheckAssociation';
switch(sim)
    case 'InterferersDistribution'
        params.la_s =  [1e-2 5e-2] ;
        params.la_u =   300e-6 ;               % users density (users/m2)
        
        
        params.alpha = 3;
        pdf = get_the_PDF_of_Interferers(params);
    case 'CompareToRef'
        params.channel = 'Rayleigh';
        params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_u =   500e-6 ;               % users density (users/m2)
        
        
        params.alpha = 3;
        R_math  = com_downlink_rate_with_smallcell_density(params);
        R_simul = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params);
        R_ref = com_downlink_rate_with_smallcell_densityRef(params)
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,1);
        f1 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-',params.la_s,R_ref,'r:');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 3 5]),{'Simulation' , 'Corollary (2)' , 'Ref. [??]'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.alpha = 4;
        R_math  = com_downlink_rate_with_smallcell_density(params);
        R_simul = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params);
        R_ref = com_downlink_rate_with_smallcell_densityRef(params)
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,2);
        f2 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-',params.la_s,R_ref,'r:');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        Hg = legend(f2([1 3 5]),{'Simulation' , 'Corollary (2)' , 'Ref. [??]'},'FontSize',25,'FontWeight','bold','Location','northwest');
        grid on;
    case 'CheckAssociation'
        params.M = 2;
        params.la_s = 1e-3;
        params.la_u =   200e-6 ;
        checkAssociation(params , index)
    case 'RangeOfDistanceToNthNeighbour'
        %params.la_s =  [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1 2e-1:1e-1:1] ;
        %params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.M = 5;
        params.la_s =  [1e-3 1e-2 1e-1];
        params.la_u =   600e-6 ;               % users density (users/m2)
        params.alpha = 4;
        params.channel = 'Rayleigh';
        range = gen_range_of_distance_with_smallcell_density_ray_multiple_assoc(params);
        
    case 'OvelappedMultiCells'
        params.la_s =  [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1 2e-1:1e-1:1] ;
        %params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        
        %params.la_s =  [1e-3:1e-3:1e-2 2e-2:1e-2:1e-1] ;
        
        params.alpha = 4;
        params.channel = 'Rayleigh';
        
        params.la_u =   100e-6 ;               % users density (users/m2)
        M = 2:params.M_agg;
        
        P_overlappeng = zeros(numel(M),numel( params.la_s));
        
        for i = 1:numel(M)
            params.M = M(i);
            P_overlappeng(i,:) = gen_Pover_with_smallcell_density_ray_multiple_assoc(params);
        end
        subplot(1,2,1);
        f1 = semilogx(params.la_s ,P_overlappeng,'k-');
        title('$$\alpha = 4,~\lambda_u = 100 $$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Probability of overlapped multicells','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        
        params.la_u =   300e-6 ;               % users density (users/m2)
        
        for i = 1:numel(M)
            params.M = M(i);
            P_overlappeng(i,:) = gen_Pover_with_smallcell_density_ray_multiple_assoc(params);
        end
        subplot(1,2,2);
        f1 = semilogx(params.la_s ,P_overlappeng,'k-');
        title('$$\alpha = 4,~\lambda_u = 300$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        %legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Probability of overlapped multicells','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
    case 'CellsDensityADR'
        params.channel = 'Rayleigh';
        params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_u =   500e-6 ;               % users density (users/m2)
        
        
        params.alpha = 3;
        R_math  = com_downlink_rate_with_smallcell_density(params);
        R_simul = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params);
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,1);
        f1 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.alpha = 4;
        R_math  = com_downlink_rate_with_smallcell_density(params);
        R_simul = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params);
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,2);
        f2 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        Hg = legend(f2([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        grid on;
        strings = {'BS1' , 'BS2' , 'BS3' , 'BS4' , 'BS5'};
        %hL = legend(f2(params.M + 1:2*params.M),strings(1:params.M),'FontSize',25,'FontWeight','bold','Location','southwest');
        
    case 'CellsDensityASE'
        params.channel = 'Rayleigh';
        params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_u =   500e-6 ;               % users density (users/m2)
        
        
        params.alpha = 3;
        R_math  = com_ASE_with_smallcell_density(params);
        R_simul = gen_ASE_with_smallcell_density_ray_multiple_assoc(params);
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,1);
        f1 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        params.alpha = 4;
        R_math  = com_ASE_with_smallcell_density(params);
        R_simul = gen_ASE_with_smallcell_density_ray_multiple_assoc(params);
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,2);
        f2 = semilogx(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        Hg = legend(f2([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        grid on;
        strings = {'BS1' , 'BS2' , 'BS3' , 'BS4' , 'BS5'};
        %hL = legend(f2(params.M + 1:2*params.M),strings(1:params.M),'FontSize',25,'FontWeight','bold','Location','southwest');
        
    case 'CellsDensityASE_Aggregate'
        params.channel = 'Rayleigh';
        params.la_s =  [1e-3 3e-3 5e-3 8e-3 1e-2 3e-2 5e-2 8e-2 1e-1] ;
        params.la_u =   500e-6 ;               % users density (users/m2)
        M = 1:params.M_agg;
        
        params.alpha = 3;                      % pass loss exponent
        R_math = zeros(numel(M),numel( params.la_s));
        R_simul = zeros(numel(M),numel( params.la_s));
        for i = 1:numel(M)
            params.M = M(i);
            R_math(i,:)  = com_ASE_Aggregate_with_smallcell_density(params);
            R_simul(i,:) = gen_ASE_Aggregate_with_smallcell_density_ray_multiple_assoc(params);
        end
        
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,1);
        f1 = semilogx(params.la_s ,R_simul/1e3,'ko' ,params.la_s,R_math/1e3,'k-');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (knats/s/Hz/$km^2$)','Interpreter','LaTex');
        
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        
        
        params.alpha = 4;                      % pass loss exponent
        R_math = zeros(numel(M),numel( params.la_s));
        R_simul = zeros(numel(M),numel( params.la_s));
        for i = 1:numel(M)
            params.M = M(i);
            R_math(i,:)  = com_ASE_Aggregate_with_smallcell_density(params);
            R_simul(i,:) = gen_ASE_Aggregate_with_smallcell_density_ray_multiple_assoc(params);
        end
        %R_math_ref  = com_downlink_rate_with_smallcell_densityRef(params);
        
        
        subplot(1,2,2);
        f2 = semilogx(params.la_s ,R_simul/1e3,'ko' ,params.la_s,R_math/1e3,'k-');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('Small cell density (cells/$m^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (knats/s/Hz/$km^2$)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        Hg = legend(f2([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        grid on;
        ax.YAxis.Exponent = 3;
        strings = {'BS1' , 'BS2' , 'BS3' , 'BS4' , 'BS5'};
        %hL = legend(f2(params.M + 1:2*params.M),strings(1:params.M),'FontSize',25,'FontWeight','bold','Location','southwest');
    case 'UsersDensityASE'
        params.channel = 'Rayleigh';
        params.la_s =  1e-2 ;
        params.la_u =   (100:100:1000) * 1e-6 ;               % users density (users/m2)
        
        
        params.alpha = 3;                      % pass loss exponent
        
        R_math  = com_ASE_with_user_density(params);
        R_simul = gen_ASE_with_user_density_ray_multiple_assoc(params);
        
        
        
        subplot(1,2,1);
        f1 = plot(params.la_u*1e6 ,R_simul,'ko' ,params.la_u*1e6 ,R_math,'k-');
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('User density (users/$km^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        xlim([100 1000]);
        
        params.alpha = 4;
        R_math  = com_ASE_with_user_density(params);
        R_simul = gen_ASE_with_user_density_ray_multiple_assoc(params);
        
        subplot(1,2,2);
        f2 = plot(params.la_u*1e6 ,R_simul,'ko' ,params.la_u*1e6,R_math,'k-');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('User density (users/$km^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (nats/s/Hz/$km^2$)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        legend(f2([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlim([100 1000]);
        grid on;
        
    case 'UsersDensityADR'
        params.channel = 'Rayleigh';
        params.la_s =  1e-2 ;
        params.la_u =   (100:100:1000) * 1e-6 ;               % users density (users/m2)
        
        
        params.alpha = 3;                      % pass loss exponent
        
        R_math  = com_downlink_rate_with_user_density(params);
        R_simul = gen_downlink_rate_with_user_density_ray_multiple_assoc(params);
        
        
        
        subplot(1,2,1);
        f1 = plot(params.la_u*1e6 ,R_simul,'ko' ,params.la_u*1e6 ,R_math,'k-');
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northeast');
        xlabel('User density (users/$km^2$)  ','Interpreter','LaTex');
        ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        xlim([100 1000]);
        
        params.alpha = 4;
        R_math  = com_downlink_rate_with_user_density(params);
        R_simul = gen_downlink_rate_with_user_density_ray_multiple_assoc(params);
        
        subplot(1,2,2);
        f2 = plot(params.la_u*1e6 ,R_simul,'ko' ,params.la_u*1e6,R_math,'k-');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('User density (cells/$km^2$)  ','Interpreter','LaTex');
        ylabel('Average downlink rate (nats/s/Hz)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        legend(f2([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northeast');
        xlim([100 1000]);
        grid on;
        
    case 'UsersDensityASE_Aggregate'
        params.channel = 'Rayleigh';
        params.la_s =  1e-2 ;
        params.la_u =   (100:100:1000) * 1e-6 ;               % users density (users/m2)
        M = 1:params.M_agg;
        
        params.alpha = 3;                      % pass loss exponent
        R_math = zeros(numel(M),numel( params.la_u));
        R_simul = zeros(numel(M),numel( params.la_u));
        for i = 1:numel(M)
            params.M = M(i);
            R_math(i,:)  = com_ASE_Aggregate_with_user_density(params);
            R_simul(i,:) = gen_ASE_Aggregate_with_user_density_ray_multiple_assoc(params);
        end
        
        
        subplot(1,2,1);
        f1 = plot(params.la_u*1e6 ,R_simul/1e3,'ko' ,params.la_u*1e6 ,R_math/1e3,'k-');
        set(f1,'MarkerSize',15);
        set(f1,'LineWidth',4);
        legend(f1([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlabel('User density (users/$km^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (knats/s/Hz/$km^2$)','Interpreter','LaTex');
        title('$$\alpha = 3$$','Interpreter','LaTex')
        grid on;
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        xlim([100 1000]);
        
        params.alpha = 4;
        for i = 1:numel(M)
            params.M = M(i);
            R_math(i,:)  = com_ASE_Aggregate_with_user_density(params);
            R_simul(i,:) = gen_ASE_Aggregate_with_user_density_ray_multiple_assoc(params);
        end
        
        subplot(1,2,2);
        f2 = plot(params.la_u*1e6 ,R_simul/1e3,'ko' ,params.la_u*1e6,R_math/1e3,'k-');
        title('$$\alpha = 4$$','Interpreter','LaTex')
        set(f2,'MarkerSize',15);
        set(f2,'LineWidth',4);
        xlabel('User density (users/$km^2$)  ','Interpreter','LaTex');
        ylabel('Area spectral efficiency (knats/s/Hz/$km^2$)','Interpreter','LaTex');
        set(gca, 'FontSize', 30);
        set(gca, 'FontWeight', 'Bold');
        legend(f2([1 params.M+1]),{'Simulation' , 'Corollary (2)'},'FontSize',25,'FontWeight','bold','Location','northwest');
        xlim([100 1000]);
        grid on;
    case 'AggregateRate'
    case 'ImpactOfFading'
        
end
end

%% PDF of the interfering BSs point process

function pdf = get_the_PDF_of_Interferers(params)
N_regions =8;  % in any dimension
lin = linspace(params.simulation_area_side(1), params.simulation_area_side(2), N_regions+1);
[X,Y] = meshgrid(lin);
Y = flipud(Y);

% Divide the area
n = 0;
rectX  = zeros(N_regions^2 , 5);
rectY  = zeros(N_regions^2 , 5);
for i = 1:N_regions
    for j = 1:N_regions
        n = n + 1;
        Xsub = X(i:i+1,i:i+1);
        Xr = reshape(Xsub,[1 4]);
        Xo = Xr([1 2 4 3]);
        Xcir = [Xo Xo(1)];
        rectX(n,:) =  Xcir;
        
        Ysub = Y(j:j+1,j:j+1);
        Yr = reshape(Ysub,[1 4]);
        Yo = Yr([1 2 4 3]);
        Ycir = [Yo Yo(1)];
        rectY(n,:) =  Ycir;    
    end
    
    
end

% for n = 1:N_regions^2
% plot(rectX(n,:),rectY(n,:),'LineWidth',2) % polygon
% end



simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
N_points = zeros(points, N_regions^2 * params.space_realizations);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        %[r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        % distribution of interferers on the regions
        N_inter = zeros(1,N_regions^2);
        for i = 1:N_users;
            server = Server_u(i,:);
            Interferers =  Server_u(Server_u ~= server); 
            for k = 1:numel(Interferers)
                xi = cells_pos(Interferers(k),1);
                yi = cells_pos(Interferers(k),2);
                
                %plot(xi,yi,'ro')
                for r = 1 : N_regions^2
                    [in,on] = inpolygon(xi,yi,rectX(r,:),rectY(r,:));
                      if (in == 1 || on == 1) 
                          N_inter(r) = N_inter(r) + in + on;
                          break;
                      end
                    
                end
            end
            
        end      
        
        
        
        
     N_points(p,N_regions^2 * (m-1) +1: m*N_regions^2)  =  N_inter ; 
    end
    
end

save  inter_pdf.mat params N_points;

pdf = N_points;
end

%% Range of Distance to Nth Neighbour
function [range_min range_max] = gen_range_of_distance_with_smallcell_density_ray_multiple_assoc(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Range_U_min = zeros(points,params.M);
Range_U_max = zeros(points,params.M);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    
    mini = 1000;
    maxi = 0;
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Distances_Users(i,:) = r_servers;
        end
        min_ = min(Distances_Users);
        max_ = max(Distances_Users);
        if(min_ < mini)
            mini = min_;
        end
        if(max_ > maxi)
            maxi = max_;
        end
    end
    Range_U_min(p,:) = mini ;
    Range_U_max(p,:) = maxi ;
    
end
range_min = Range_U_min;
range_max = Range_U_max;
save rangef.mat params Range_U_min Range_U_max ;
end
%% Probability of overlapping
function [P] = gen_Pover_with_smallcell_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
P_U = zeros(points , params.space_realizations);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      Network Plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  s = plot(cells_pos(:,1),cells_pos(:,2),'sb');
        %  for i = 1:N_cells
        %     text(cells_pos(i,1) - 2,cells_pos(i,2),num2str(i),'Color','blue','FontSize',14);
        %  end
        %  set(s,'MarkerSize',20);
        %  xlabel('x (meters)');
        %  ylabel('y (meters)');
        %  set(gca, 'FontName', 'Arial');
        %  set(gca, 'FontSize', 20);
        %  set(gca, 'FontWeight', 'Bold');
        %
        %  hold on;
        %
        %  u = plot(users_pos(:,1),users_pos(:,2),'or');
        %  set(u,'MarkerSize',15);
        %  for i = 1:N_users
        %      text(users_pos(i,1) - 2,users_pos(i,2),num2str(i),'Color','red','FontSize',14);
        %  end
        %
        %
        %
        %
        %  for i = 1:N_users;
        %      for s = 1:params.M
        %          x = [users_pos(i,1) , cells_pos(Server_u(i,s),1) ];
        %          y = [users_pos(i,2) , cells_pos(Server_u(i,s),2) ];
        %
        %          line(x ,y,'Color','red','LineStyle','-','LineWidth',3)
        %      end
        %
        %  end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this is an approximate indicator, we have to check it in another
        % way
        overlapping_indicator = zeros(N_cells,1);
        Aconnections = sum(A');   % if a cell serves more than user the correspondong element will be more than 1
        overlapping_indicator = (Aconnections > 1);
        Active = (Aconnections >= 1);
        N_ovelapped_multicells = sum (overlapping_indicator);
        N_active_cells = sum(Active);
        P_U(p,m) = N_ovelapped_multicells / N_active_cells;
        % detect ovelapped multicells
        %        overlapping_indicator = zeros(N_users,1);
        %         for i = 1:N_users;
        %             for s = 1:params.M
        %                 server = Server_u(i,s);
        %                 if (sum(A(server,:)) > 1)
        %                     ind = find(A(server,:) == 1);
        %                     overlapping_indicator(ind) = 1;
        %                 end
        %             end
        %
        %         end
        %
        %P_U(p,m) = sum(overlapping_indicator) / N_users;
    end
    
end
normfact = params.space_realizations;


P = sum(P_U,2) / normfact;
end

%% Average Downlink Rate (ADR) with User Density
function [R_math] = com_downlink_rate_with_user_density(params)

points = numel(params.la_u);

R_math = zeros(points,params.M);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    m = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    k = params.la_s / params.la_u(p);   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
    factor = 1;
    switch(params.M)
        case 1
            F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 2
            F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 3
            F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
            
        case 4
            F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
        otherwise
            disp('Invalid Multicell size (M)')
    end
    
end

end
function [R] = gen_downlink_rate_with_user_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_u);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Smallcell Density: ' , num2str(params.la_s)]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Users Density: ' , num2str(params.la_u(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u(p) * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = (Int_all - Int_servers) / params.M; % this is modified (remove  params.M)
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(Rate_U,4),3) / normfact;
end

%% Area Spectral Efficiency (ASE) with User Density
function [R_math] = com_ASE_with_user_density(params)

points = numel(params.la_u);

R_math = zeros(points,params.M);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    m = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    k = params.la_s / params.la_u(p);   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
    factor = 1;
    switch(params.M)
        case 1
            F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 2
            F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 3
            F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
            
        case 4
            F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
        otherwise
            disp('Invalid Multicell size (M)')
    end
    
end

end
function [R] = gen_ASE_with_user_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_u);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Smallcell Density: ' , num2str(params.la_s)]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Users Density: ' , num2str(params.la_u(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u(p) * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = (Int_all - Int_servers) / params.M; % this is modified (remove  params.M)
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(Rate_U,4),3) / normfact;
end

%% Area Spectral Efficiency (ASE) with User Density
% function [R_math] = com_ASE_with_user_density(params)
%
% points = numel(params.la_u);
%
% R_math = zeros(points,params.M);
%
% % Rice Channel Parameters
% if strcmp(params.channel ,'Rice')
%     K = params.rice.K;
%     m = (K+1)^2/(2*K+1);
%     sigma = 1;
%     omega = (2*K+1) * sigma^2;
% end
%
% for p = 1:points
%     k = params.la_s / params.la_u(p);   % densification ratio
%     po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
%
%     switch(params.M)
%         case 1
%             F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
%         case 2
%             F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%
%             R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%             F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
%         case 3
%             F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%
%             R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%             F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%             F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%         case 4
%             F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%
%             R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%             F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%             F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
%
%             F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
%                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
%             R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
%         otherwise
%             disp('Invalid Multicell size (M)')
%     end
%     simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
%     R_math(p,:) = R_math(p,:) .* params.la_u(p) * simulation_area  ;
% end
%
% end
% function [R] = gen_ASE_with_user_density_ray_multiple_assoc(params)
%
% simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
% points = numel(params.la_u);
% Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
% for p = 1:points
%     fprintf('\n')
%     disp(['Multicell Size: ' , num2str(params.M)]);
%     disp(['Fading Channel: ' , params.channel]);
%     disp(['Smallcell Density: ' , num2str(params.la_s)]);
%     disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
%     disp(['Users Density: ' , num2str(params.la_u(p))]);
%     disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
%     disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
%     for m = 1:params.space_realizations;
%         if(mod(m,params.space_realizations/100) == 0)
%             fprintf('|');
%         end
%         mu_s = params.la_s * simulation_area;
%         mu_u = params.la_u(p) * simulation_area;
%
%         N_cells = poissrnd(mu_s);
%         N_users = poissrnd(mu_u);
%
%         cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
%         users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
%
%         r_su = pdist2(cells_pos,users_pos,'euclidean') ;
%
%         A = zeros(N_cells,N_users);  % association matrix
%         S_u = zeros(N_users,params.M);
%         I_u = zeros(N_users,params.M);
%         SIR_u = zeros(N_users,params.M);
%         Server_u = zeros(N_users,params.M);
%         r_u = zeros(N_users,params.M);
%         % Association  [Long_Term Association i.e. the channel is averaged to 1]
%         [r_u_ , Server_u_] = min(r_su);
%         for i = 1:N_users;
%             [servers, r_servers] = getServers(i,r_su, params.M);
%             Server_u(i,:) = servers;
%             r_u(i,:) = r_servers;
%             A(servers,i) = 1;
%         end
%
%         for t = 1:params.time_slots
%             Hui = exprnd(1,N_cells,N_users) ;
%             R_u =  params.P *  Hui.* r_su.^-params.alpha;
%             Ac = sum(A,2)>=1;
%
%             for i = 1:N_users;
%                 for s = 1:params.M
%                     S_u(i,s) = R_u(Server_u(i,s),i);
%                 end
%             end
%             Int_all = sum(bsxfun(@times,Ac,R_u));
%             Int_servers = sum(S_u,2)';
%             I_u = Int_all - Int_servers;
%
%             for i = 1:N_users;
%                 for s = 1:params.M
%                     SIR_u(i,s) = S_u(i,s) / I_u(i);
%                 end
%             end
%             Rate_u = log(1 + SIR_u);
%             %Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
%             Rate_U(p,:,m,t) = sum(Rate_u);
%         end
%     end
%
% end
% normfact = params.space_realizations * params.time_slots ;% * simulation_area;
% Rate_U(isinf(Rate_U)) = 0;
% Rate_U(Rate_U == 0) = max(Rate_U(:));
% R = sum(sum(Rate_U,4),3) / normfact;
% end

%% Aggregate Area Spectral Efficiency with User Density
function [R_math_aggregate] = com_ASE_Aggregate_with_user_density(params)

points = numel(params.la_u);

R_math_aggregate = zeros(1,points);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    m = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    k = params.la_s / params.la_u(p);   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
    factor = 1;
    F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo  + (params.M-1) .* hyp2f1(1,2/a,1+ 2/a,-1./z ))./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
        mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    R_math_aggregate(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
    %     switch(params.M)
    %         case 1
    %             F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
    %         case 2
    %             F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %
    %             R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
    %         case 3
    %             F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %
    %             R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %         case 4
    %             F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %
    %             R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
    %             R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
    %         otherwise
    %             disp('Invalid Multicell size (M)')
    %     end
    simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    %     R_math(p,:) = R_math(p,:) .* params.la_u(p) * simulation_area  ;
    %     R_math_aggregate(p) = sum(R_math(p,:));
    R_math_aggregate(p) = R_math_aggregate(p) .* params.la_u(p) * 1e6 ; %/ params.M;
end

end

function [R] = gen_ASE_Aggregate_with_user_density_ray_multiple_assoc_old(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_u);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Smallcell Density: ' , num2str(params.la_s)]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Users Density: ' , num2str(params.la_u(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u(p) * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = Int_all - Int_servers;
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            %Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_U(p,:,m,t) = sum(Rate_u);
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;
Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(sum(Rate_U,4),3),2) / normfact;
end

function [R] = gen_ASE_Aggregate_with_user_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_u);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Users Density: ' , num2str(params.la_u(p))]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Smallcell Density: ' , num2str(params.la_s)]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s * simulation_area;
        mu_u = params.la_u(p) * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = (Int_all - Int_servers) / params.M; % this is modified (remove  params.M)
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            %Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_U(p,:,m,t) = sum(Rate_u);
        end
    end
    
end
normfact = params.space_realizations * params.time_slots *( simulation_area * 1e-6) ; %* params.M ;

Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(sum(Rate_U,4),3),2) / normfact;
end

%% Average Downlink Rate (ADR) with Small cell Density
function [R_math] = com_downlink_rate_with_smallcell_density(params)

points = numel(params.la_s);

R_math = zeros(points,params.M);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    mo = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    %factor = params.M;
    factor = 1;
    switch(params.M)
        case 1
            F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 2
            F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 3
            F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
            
        case 4
            F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
        otherwise
            disp('Invalid Multicell size (M)')
    end
    
end

end
function [R] = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = (Int_all - Int_servers) / params.M; % this is modified (remove  params.M)
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;

Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(Rate_U,4),3) / normfact;
end

%% Area Spectral Efficiency (ASE) with Small cell Density
function [R_math] = com_ASE_with_smallcell_density(params)

points = numel(params.la_s);

R_math = zeros(points,params.M);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    mo = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    %factor = params.M;
    factor = 1;
    switch(params.M)
        case 1
            F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 2
            F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
        case 3
            F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
            
        case 4
            F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
            
            F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
                mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
            R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
        otherwise
            disp('Invalid Multicell size (M)')
    end
    simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    R_math(p,:) = R_math(p,:) .* params.la_u * simulation_area  ;
    
end

end
function [R] = gen_ASE_with_smallcell_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = (Int_all - Int_servers) / params.M; % this is modified (remove  params.M)
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            %Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_U(p,:,m,t) = sum(Rate_u);
        end
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;

Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(Rate_U,4),3) / normfact;
end

%% Aggregate Area Spectral Efficiency (ASE) with Small cell Density
function [R_math_aggregate] = com_ASE_Aggregate_with_smallcell_density(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);

R_math_aggregate = zeros(1,points);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    mo = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    %factor = params.M;
    factor = 1;
    F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo  + (params.M-1) .* hyp2f1(1,2/a,1+ 2/a,-1./z ))./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
        mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    R_math_aggregate(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
    %     switch(params.M)
    %         case 1
    %             F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf);
    %         case 2
    %             F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %
    %             R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf);
    %         case 3
    %             F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %
    %             R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %         case 4
    %             F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %
    %             R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf);
    %
    %             F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + factor.*z./mi).^-mi + ...
    %                 mi.*(mi).^mi.*(1 - 2/a).^-1 .* factor.* z .* (factor.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,factor.*z.*(factor.*z+mi).^-1))).^n) ;
    %             R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf);
    %         otherwise
    %             disp('Invalid Multicell size (M)')
    %     end
    %simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
    %R_math(p,:) = R_math(p,:) .* params.la_u * simulation_area  ;
    %R_math_aggregate(p) = sum(R_math(p,:));
    R_math_aggregate(p) = R_math_aggregate(p) .* params.la_u * 1e6 ; %/ params.M;
end

end
function [R] = gen_ASE_Aggregate_with_smallcell_density_ray_multiple_assoc(params)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;
points = numel(params.la_s);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
for p = 1:points
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
    disp(['Users Density: ' , num2str(params.la_u)]);
    disp(['Fading Channel: ' , params.channel]);
    disp(['Path Loss Exponent: ' , num2str(params.alpha)]);
    disp(['Smallcell Density: ' , num2str(params.la_s(p))]);
    disp('          10%       20%       30%       40%       50%       60%       70%       80%       90%       100%');
    disp('|.........|.........|.........|.........|.........|.........|.........|.........|.........|.........|');
    for m = 1:params.space_realizations;
        if(mod(m,params.space_realizations/100) == 0)
            fprintf('|');
        end
        mu_s = params.la_s(p) * simulation_area;
        mu_u = params.la_u * simulation_area;
        
        N_cells = poissrnd(mu_s);
        N_users = poissrnd(mu_u);
        
        cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
        users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);
        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
        
        A = zeros(N_cells,N_users);  % association matrix
        S_u = zeros(N_users,params.M);
        I_u = zeros(N_users,params.M);
        SIR_u = zeros(N_users,params.M);
        Server_u = zeros(N_users,params.M);
        r_u = zeros(N_users,params.M);
        % Association  [Long_Term Association i.e. the channel is averaged to 1]
        [r_u_ , Server_u_] = min(r_su);
        for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            r_u(i,:) = r_servers;
            A(servers,i) = 1;
        end
        
        for t = 1:params.time_slots
            Hui = exprnd(1,N_cells,N_users) ;
            R_u =  params.P *  Hui.* r_su.^-params.alpha;
            Ac = sum(A,2)>=1;
            
            for i = 1:N_users;
                for s = 1:params.M
                    S_u(i,s) = R_u(Server_u(i,s),i);
                end
            end
            Int_all = sum(bsxfun(@times,Ac,R_u));
            Int_servers = sum(S_u,2)';
            I_u = (Int_all - Int_servers) / params.M; % this is modified (remove  params.M)
            
            for i = 1:N_users;
                for s = 1:params.M
                    SIR_u(i,s) = S_u(i,s) / I_u(i);
                end
            end
            Rate_u = log(1 + SIR_u);
            %Rate_U(p,:,m,t) = sum(Rate_u) / N_users;
            Rate_U(p,:,m,t) = sum(Rate_u);
        end
    end
    
end
normfact = params.space_realizations * params.time_slots *( simulation_area * 1e-6);% *  params.M ;

Rate_U(isinf(Rate_U)) = 0;
Rate_U(Rate_U == 0) = max(Rate_U(:));
R = sum(sum(sum(Rate_U,4),3),2) / normfact;
end
%%
function [R_math_ref] = com_downlink_rate_with_smallcell_densityRef(params)

points = numel(params.la_s);

R_math_ref = zeros(points,1);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K;
    mo = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end

for p = 1:points
    
    
    F = @(z,a,mo,omega,mi,n) (1 - (1 + z.*omega./mo).^-mo)./(z .* ((1 + z./mi).^-mi + ...
        mi.*(mi).^mi.*(1 - 2/a).^-1 .*  z .* (z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1)).^n) ;
    R_math_ref(p) = integral(@(z)F(z,params.alpha,1,1,1,1),0,inf);
    
    
end

end

function [servers,r_servers] = getServers(user,distances, N)
servers = zeros(N,1);
r_servers = zeros(N,1);
for i = 1:N
    [r_servers(i) , servers(i)] = min(distances(:,user));
    distances(servers(i),user) = inf;
end
end

function checkAssociation(params,index)

simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;

mu_s = params.la_s(index) * simulation_area;
mu_u = params.la_u * simulation_area;


N_cells = poissrnd(mu_s);
N_users = poissrnd(mu_u);

cmap = hsv(N_users);

cells_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
users_pos = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users, 2);


r_su = pdist2(cells_pos,users_pos,'euclidean') ;


A = zeros(N_cells,N_users);  % association matrix


for i = 1:N_users;
    [servers, r_servers] = getServers(i,r_su, params.M);
    Server_u(i,:) = servers;
    r_u(i,:) = r_servers;
    A(servers,i) = 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Network Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = plot(cells_pos(:,1),cells_pos(:,2),'sb');
for i = 1:N_cells
    text(cells_pos(i,1) - 2,cells_pos(i,2),num2str(i),'Color','blue','FontSize',14);
end
set(s,'MarkerSize',20);
xlabel('x (meters)');
ylabel('y (meters)');
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 20);
set(gca, 'FontWeight', 'Bold');

hold on;

u = plot(users_pos(:,1),users_pos(:,2),'or');
set(u,'MarkerSize',15);
for i = 1:N_users
    text(users_pos(i,1) - 2,users_pos(i,2),num2str(i),'Color','red','FontSize',14);
end




for i = 1:N_users;
    for s = 1:params.M
        x = [users_pos(i,1) , cells_pos(Server_u(i,s),1) ];
        y = [users_pos(i,2) , cells_pos(Server_u(i,s),2) ];
        
        line(x ,y,'Color','red','LineStyle','-','LineWidth',3)
    end
    
end
end