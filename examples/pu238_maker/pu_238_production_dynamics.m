%% Simple Estimate of Np-237 -> Pu-238
clear
clc
close 'all'

%%
% 
%  The purpose of this script is to create a simple and yet quantitatively
%  representative estimate of the production/decay chains of isotopes
%  relevant to Pu-238 production from Np-237 targets.  
% 

%% Parameters

flux = 1e14; % n/cm^2-s, assumed constant
time = 1000*3600; % s, 1000 hr converted to seconds
num_ts = 1e6; % number of time steps in the discrete simulation
tSpace = linspace(0,time,num_ts);


N237 = 0.004988; % atom/b-cm, atom density of Np-237 in target.
% estimated from NpO2/Aluminum CERMET target properties
N237 = N237*1e24; % conver to atoms/cm^3

%% Isotopes to track:
% 1 - Np-237
% 2 - Np-238
% 3 - Pu-238
% 4 - Pu-239
% 5 - Pu-240
% 6 - Pu-241
% 7 - Np-239 <-- added later
N_isotopes = 7;
%% Nuclear Data
sigma_c = nan(N_isotopes,1);% barns, n,gamma thermal cross sections 
sigma_f = nan(N_isotopes,1);% barns, fission thermal cross sections
lam = nan(N_isotopes,1); % sec, decay constants

% 1 - Np-237
sigma_c(1) = 169;
sigma_f(1) = 0;
lam(1) = 0; % assume stable

% 2 - Np-238
sigma_c(2) = 480;
sigma_f(2) = 2100;
lam(2) = log(2)/(50.88*3600); 

% 3 - Pu-238
sigma_c(3) = 540;
sigma_f(3) = 18; % use data from chart of nuclides
lam(3) = 0;

% 4 - Pu-239
sigma_c(4) = 271;
sigma_f(4) = 750;
lam(4) = 0;

% 5 - Pu-240
sigma_c(5) = 290;
sigma_f(5) = 0;
lam(5) = 0;

% 6 - Pu-241
sigma_c(6) = 361;
sigma_f(6) = 1010;
lam(6) = 0;

% 7 - Np-239 (added last; sorry)
sigma_c(7) = 60;
sigma_f(7) = 0;
lam(7) = log(2)/(2.35*24*3600);

%% Form Matrix for Bateman Equation
M = zeros(N_isotopes,N_isotopes);
M(1,1) = -(sigma_c(1) + sigma_f(1));
M(2,2) = -(sigma_c(2) + sigma_f(2)); M(2,1) = sigma_c(1);
M(3,3) = -(sigma_c(3) + sigma_f(3)); 
M(4,4) = -(sigma_c(4) + sigma_f(4)); M(4,3) = sigma_c(3);
M(5,5) = -(sigma_c(5) + sigma_f(5)); M(5,4) = sigma_c(4);
M(6,6) = -(sigma_c(6) + sigma_f(6)); M(6,5) = sigma_c(5);
M(7,7) = -(sigma_c(7) + sigma_f(7)); M(7,2) = sigma_c(2);
M = M*(1e-24); % convert to cm^2

lam = diag(lam); % convert to matrix with lam along diagonal.
lam(3,2) = -lam(2,2); % ugly but correctly expresses Pu-238 prod from Np-238 decay.
lam(4,7) = -lam(7,7); % again, ugly but gives Pu-239 prod from decay of Np-239
%% Pre-allocate Data Array
N = zeros(N_isotopes,num_ts); % store isotope inventory over time
N(1,1) = N237; % initialize with starting atom density of Np-237
dt = time/(num_ts - 1); % sec, time step size

%% Commence time stepping

for t = 1:(num_ts - 1)

    if mod(t,100000) == 0
        fprintf('Commencing time step %g.\n',t);
    end
    dN = ((M*N(:,t))*flux - lam*N(:,t))*dt;
    N(:,t+1) = N(:,t) + dN;
end

% get relative atom fractions
N_rel = N/N(1,1); 

% estimate Pu-238 21 days after exposure
% assume all Np-238 decays in that time fram.
N_p21 = N(3,:) + N(2,:);
N_p21_rel = N_p21/N(1,1);

%% Plot some results
figure(1)
tSpace = tSpace/3600; % convert to hours.
mark_ind = 1:100000:num_ts;

title('System Dynamics of Pu-238 Production',...
    'fontsize',16,'fontweight','bold');
xlabel('Irradiation time [hr]','fontsize',14,'FontWeight','bold');


yyaxis left
semilogy(tSpace,N_rel(1,:),'-ro','linewidth',2,...
    'MarkerIndices',mark_ind);
hold on;
semilogy(tSpace,N_rel(2,:),'-bs','linewidth',2,...
    'MarkerIndices',mark_ind);
semilogy(tSpace,N_rel(3,:),'-gd','linewidth',2,...
    'MarkerIndices',mark_ind);
semilogy(tSpace,N_rel(4,:),'--ms','linewidth',2,...
    'MarkerIndices',mark_ind);
semilogy(tSpace,N_rel(7,:),'-ko','linewidth',2,...
    'MarkerIndices',mark_ind);
semilogy(tSpace,N_p21_rel,'-.r^','linewidth',2,...
    'MarkerIndices',mark_ind);

hold off;


ylabel('Atom Density [atom/b-cm]',...
    'FontSize',14,'FontWeight','bold');
axis([0 1000 1e-6 1]);

yyaxis right

semilogy(tSpace,N_rel(3,:)./N_rel(4,:),'-gs','linewidth',2,...
    'MarkerIndices',mark_ind);
ylabel('Pu-238/Pu-239','Fontsize',14,'FontWeight','bold');



grid on
legend('Np-237','Np-238','Pu-238','Pu-239','Np-239',...
    'Pu-238 (after 21 days)','Pu-238/Pu-239');