%% U-238 cross section plot script

clear
clc
close 'all'

addpath('/opt/mcnp/endf71_h5/mcnp_endfb71');

% to see details on the contents of a H5 file, use h5disp()
%h5disp('U238.h5');

% for each nuclide there are data sets organized by temperature
% and reaction type.

% Reaction types are numbered in accordance with the ENDF MT numbers.
% Listings of these reactions can be found in various locations, for
% example: http://serpent.vtt.fi/mediawiki/index.php/ENDF_reaction_MT%27s_and_macroscopic_reaction_numbers

% for this plot, we will use the (n,gamma) reaction (reaction 102)

nuclide = 'U238';
data_file = strcat(nuclide,'.h5');
reaction = 'reaction_102';
temperatures = {'0K','250K','294K','600K','900K','1200K','2500K'};
lineSpecs = {'-y','-m','-c','-r','-g','-b','-k'};


energy_stb = strcat('/',nuclide,'/energy');
reaction_stb = strcat('/',nuclide,'/reactions/',reaction);

min_en = 6; % eV
max_en = 7.5; % eV

%for energy = 1:length(temperatures)
for energy = 1:length(temperatures)
    en = h5read(data_file,...
        strcat(energy_stb,'/',temperatures{energy}));
    
    idx = find(((en >= min_en) & (en <= max_en)));
    cross_section = h5read(data_file,...
        strcat(reaction_stb,'/',temperatures{energy},'/xs'));
    
    
    loglog(en(idx),cross_section(idx),...
        lineSpecs{energy},'linewidth',3);
    
    if energy == 1
        hold on
        
    end
    
end

title('U-238 (n,\gamma) Doppler Broadening','fontsize',18,'fontweight','bold');
xlabel('Energy (eV)','fontsize',16,'fontweight','bold');
ylabel('U-238 (n,\gamma) Cross Section (barns)','fontsize',16,'fontweight','bold');
legend(temperatures,'fontsize',16);
grid on
set(gca,'fontsize',14,'fontweight','bold');
