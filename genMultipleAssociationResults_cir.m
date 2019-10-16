function  genMultipleAssociationResults()
clc;close all; clear;
%format long e;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.channel = 'Rayleigh';
params.la_s =  [1e-3  5e-3 1e-2];% 5e-2 1e-1 5e-1 1] ; 
params.la_u =   500e-6 ;               % users density (users/m2)
params.alpha = 4;                      % pass loss exponent
params.M = 2;                                 
params.P = 1; 

params.R = 200 * ones(1,numel(params.la_s));
params.space_realizations = 5000;
params.time_slots = 10;

R_math  = com_downlink_rate_with_smallcell_density(params);
R_simul = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params);




    g = plot(params.la_s ,R_simul,'ko' ,params.la_s,R_math,'k-');
    legend(g([1 params.M+1]),{'Simulation' , 'Theorem (2)'},'FontSize',20,'FontWeight','bold');
    set(g,'MarkerSize',15);
    set(g,'LineWidth',4);
    xlabel('Small cell density (cells/m^2)  ');
    ylabel('Average downlink rate (nats/s/Hz)'); 
    set(gca, 'FontSize', 25);
    set(gca, 'FontWeight', 'Bold');
    
    grid on;

end
%%
function [R_math] = com_downlink_rate_with_smallcell_density(params)

points = numel(params.la_s);

R_math = zeros(points,params.M);

% Rice Channel Parameters
if strcmp(params.channel ,'Rice')
    K = params.rice.K; 
    m = (K+1)^2/(2*K+1);
    sigma = 1;
    omega = (2*K+1) * sigma^2;
end
    
for p = 1:points
    k = params.la_s(p) / params.la_u;   % densification ratio
    po = ((3.5 * k) ./ (1 + 3.5 * k)).^ 3.5;
    
    switch(params.M)
        case 1
            F = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p) = integral(@(z)F(z,params.alpha,1,1,1,po,params.M),0,inf); 
        case 2
            F21 = @(z,a,mo,omega,mi,po,n)  hyp2f1(1,2/a,1+ 2/a,-1./z) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F21(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
            F22 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F22(z,params.alpha,1,1,1,po,params.M),0,inf); 
        case 3
            F31 = @(z,a,mo,omega,mi,po,n)  ( 2 * hyp2f1(1,2/a,1+ 2/a,-1./z) -  hyp2f1(1,4/a,1+ 4/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F31(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
            F32 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,4/a,1+ 4/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F32(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
            F33 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F33(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
        case 4
            F41 = @(z,a,mo,omega,mi,po,n)  ( 3 * hyp2f1(1,2/a,1+ 2/a,-1./z) - 3 * hyp2f1(1,4/a,1+ 4/a,-1./z) +  hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            
            R_math(p,1) = integral(@(z)F41(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
            F42 = @(z,a,mo,omega,mi,po,n) ( 3* hyp2f1(1,4/a,1+ 4/a,-1./z) -  2*hyp2f1(1,6/a,1+ 6/a,-1./z) ) ./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p,2) = integral(@(z)F42(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
            F43 = @(z,a,mo,omega,mi,po,n) hyp2f1(1,6/a,1+ 6/a,-1./z)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p,3) = integral(@(z)F43(z,params.alpha,1,1,1,po,params.M),0,inf); 
            
            F44 = @(z,a,mo,omega,mi,po,n) (1 - (1 + z.*omega./mo).^-mo)./(z .*(po + (1-po) .* ((1 + params.M.*z./mi).^-mi + ...
                                          mi.*(mi).^mi.*(1 - 2/a).^-1 .* params.M.* z .* (params.M.*z+mi).^-(mi+1) .* hyp2f1(mi+1,1,2-2/a,params.M.*z.*(params.M.*z+mi).^-1))).^n) ;
            R_math(p,4) = integral(@(z)F44(z,params.alpha,1,1,1,po,params.M),0,inf); 
        otherwise
            disp('Invalid Multicell size (M)')
    end
end

end
function [R] = gen_downlink_rate_with_smallcell_density_ray_multiple_assoc(params)

points = numel(params.la_s);
Rate_U = zeros(points ,params.M, params.space_realizations , params.time_slots);
H = exprnd(1,1e6,1) ;
for p = 1:points
    simulation_area = pi * params.R(p)^2 ;
    fprintf('\n')
    disp(['Multicell Size: ' , num2str(params.M)]);
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
        
        cells_pos = generatePositions(params.R(p), N_cells);
        users_pos = generatePositions(params.R(p), N_users);
        N_users = N_users + 1;
        users_pos = [[0 0] ; users_pos];        
        r_su = pdist2(cells_pos,users_pos,'euclidean') ;
  
        A = zeros(N_cells,N_users);  % association matrix
    
        Server_u = zeros(N_users,params.M);
  
        % Association  [Long_Term Association i.e. the channel is averaged to 1] 
         for i = 1:N_users;
            [servers, r_servers] = getServers(i,r_su, params.M);
            Server_u(i,:) = servers;
            A(servers,i) = 1;
         end
         Ac = sum(A,2);
         Ac(Server_u(1,:)) =  0;
         Ac = Ac >= 1;
            
        for t = 1:params.time_slots
            Hui =H(t:N_cells+t-1) ;
            R_u =  params.P *  Hui.* r_su(:,1).^-params.alpha;
            I_u = sum(Ac.*R_u);              
            S_u(1,:) = R_u(Server_u(1,:));   
            SIR_u = S_u ./ I_u;
            Rate_u = log(1 + SIR_u);
            Rate_U(p,:,m,t) = Rate_u ;
        end
        
        
    end
    
end
normfact = params.space_realizations * params.time_slots ;% * simulation_area;

R = sum(sum(Rate_U,4),3) / normfact;
end
function [servers,r_servers] = getServers(user,distances, N)
    servers = zeros(N,1);
    r_servers = zeros(N,1);
    for i = 1:N
        [r_servers(i) , servers(i)] = min(distances(:,user));
        distances(servers(i),user) = inf;
    end
end

function [position] = generatePositions(R,N)
    r = R*sqrt(rand(N,1));
    theta = 2*pi*rand(N,1);
    x = r.*cos(theta);
    y = r.*sin(theta);
    position = [x y];
end