diab = importdata('data/BMR_1_diab_A.csv');
obes = importdata('data/BMR_2_obes_A.csv');

diab3 = abs(diab);
diab3(1:1+size(diab3,1):end) = 0;
obes3 = abs(obes);
obes3(1:1+size(obes3,1):end) = 0;

iter = 1000;

results_obes = [0 0 0 0];
results_diab = [0 0 0 0];
% Note: Added dens_net as 4th output of small_world_propensity function!
for k = 1:iter
    [SWP_obes3 delta_C_obes3 delta_L_obes3 dens_net_obes3] = small_world_propensity(obes3);
    [SWP_diab3 delta_C_diab3 delta_L_diab3 dens_net_diab3] = small_world_propensity(diab3);
    results_obes = [results_obes; [SWP_obes3 delta_C_obes3 delta_L_obes3 dens_net_obes3]];
    results_diab = [results_diab; [SWP_diab3 delta_C_diab3 delta_L_diab3 dens_net_diab3]];
end
results_obes = results_obes(2:iter+1,:);
results_diab = results_diab(2:iter+1,:);
results_obes_mean = mean(results_obes);
results_diab_mean = mean(results_diab);

csvwrite('data/SWP_diab.csv', results_diab);
csvwrite('data/SWP_obes.csv', results_obes);
