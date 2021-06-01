%% load models
clc; clear all; close all;
olddir = pwd;
cd('/cm/shared/uniol/software/IQMtools/1.2.2.2-Matlab-2017b');
installIQMtools
cd(olddir);
% set random number generator
rng('shuffle');
% Test: less effector
modelfrontend = IQMmodel('Beelen_small_2PDE_7_stoch_CSM.txtbc');

%% make frontend model irreversible
modelfrontend = IQMmakeirreversible(modelfrontend);

%% stochastic simulation frontend to get E*
V = 1e-6;
timeend = 3;
runs = 100;
Nsample = 1;
stochsimfrontend = IQMstochsim2(modelfrontend,V,timeend,runs,Nsample); 

index1 = stateindexIQM(modelfrontend,'PDE_a_Ga_GTP');
index2 = stateindexIQM(modelfrontend,'Ga_GTP_a_PDE_a_Ga_GTP');
index3 = stateindexIQM(modelfrontend, 'Ga_GTP_PDE_a_Ga_GTP');
PDE_singcell = cell(1,stochsimfrontend.runs);
PDE_doubcell = cell(1,stochsimfrontend.runs);
timefrontstochALL = cell(1,stochsimfrontend.runs);
for k=1:stochsimfrontend.runs
    PDE_sing = stochsimfrontend.speciesdata{k}(:,index1) + stochsimfrontend.speciesdata{k}(:,index3);
    PDE_doub = stochsimfrontend.speciesdata{k}(:,index2);
    timestochk = stochsimfrontend.time{k};
    PDE_singcell{k} = PDE_sing;
    PDE_doubcell{k} = PDE_doub;
    timefrontstochALL{k} = timestochk;
end
clearvars stochsimfrontend Estochk timestochk modelfrontend;