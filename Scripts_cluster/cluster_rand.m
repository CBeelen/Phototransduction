%% load models
clc; clear all; close all;
olddir = pwd;
cd('/cm/shared/uniol/software/IQMtools/1.2.2.2-Matlab-2017b');
installIQMtools
cd(olddir);
% set random number generator
rng('shuffle');
% Test: less effector
modelfrontend = IQMmodel('Beelen_small_2PDE_7_stoch.txtbc');
modelfrontend_prec = IQMmodel('Beelen_small_2PDE_7_stoch_prec.txtbc');

%% make frontend model irreversible
modelfrontend = IQMmakeirreversible(modelfrontend);
modelfrontend_prec = IQMmakeirreversible(modelfrontend_prec);

%% stochastic simulation frontend to get E*
V = 1e-6;
timeend = 3;
Nsample = 1;
runs = 100;      % number of independent simulation runs
prob = 0.018;     % probability to activate a Rh_G complex vs a single Rh
stochsimfrontend = cell(1,runs);
index1 = stateindexIQM(modelfrontend,'PDE_a_Ga_GTP');
index2 = stateindexIQM(modelfrontend,'Ga_GTP_a_PDE_a_Ga_GTP');
index3 = stateindexIQM(modelfrontend, 'Ga_GTP_PDE_a_Ga_GTP');
total = zeros(1,runs);  % total number of photoisomerizations
R0 = zeros(1,runs);     % total number of activated R
G_R0 = zeros(1,runs);   % total number of activated R*G
for i=1:runs   % for each run: num random numbers and new simulation
    total(i) = poissrnd(1);  % number of activated (total)
    single_sims = cell(1,total(i));    % cell for single simulations
    R0(i) = 0;     % activated single R0
    G_R0(i) = 0;   % activated R0_G complexes
    for j=1:total(i)
        rnd = rand;    % random number between 0 and 1 (uniformly distr.)
        if(rnd < prob) % choose R0_G compley with probability prob
            G_R0(i) = G_R0(i) + 1;
        else
            R0(i) = R0(i) + 1;    % else choose single R0
        end
    end        
    for j=1:G_R0(i)    % R is precoupled
        single_sims{j} = IQMstochsim2(modelfrontend_prec, V, timeend, 1, Nsample);
    end
    for j=1:R0(i)
        single_sims{G_R0(i)+j} = IQMstochsim2(modelfrontend, V, timeend, 1, Nsample);
    end
    for j=1:total(i)
        PDE_sing{i,j} = single_sims{j}.speciesdata{1,1}(:,index1) + single_sims{j}.speciesdata{1,1}(:,index3);
        PDE_doub{i,j} = single_sims{j}.speciesdata{1,1}(:,index2);
        time{i,j} = single_sims{j}.time{1,1};
    end
    clearvars single_sims;
end
