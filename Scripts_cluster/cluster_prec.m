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
runs = 100;
Nsample = 1;
stochsimfrontend = IQMstochsim2(modelfrontend,V,timeend,runs,Nsample); 

index1 = stateindexIQM(modelfrontend,'PDE_a_Ga_GTP');
index2 = stateindexIQM(modelfrontend,'Ga_GTP_a_PDE_a_Ga_GTP');
index3 = stateindexIQM(modelfrontend, 'Ga_GTP_PDE_a_Ga_GTP');
indexGt0 = stateindexIQM(modelfrontend,'R0_Gt');
indexGt1 = stateindexIQM(modelfrontend,'R1_Gt');
indexGt2 = stateindexIQM(modelfrontend,'R2_Gt');
indexGt3 = stateindexIQM(modelfrontend,'R3_Gt');
indexG0 = stateindexIQM(modelfrontend,'R0_G');
indexG1 = stateindexIQM(modelfrontend,'R1_G');
indexG2 = stateindexIQM(modelfrontend,'R2_G');
indexG3 = stateindexIQM(modelfrontend,'R3_G');
indexG0GTP = stateindexIQM(modelfrontend,'R0_G_GTP');
indexG1GTP = stateindexIQM(modelfrontend,'R1_G_GTP');
indexG2GTP = stateindexIQM(modelfrontend,'R2_G_GTP');
indexG3GTP = stateindexIQM(modelfrontend,'R3_G_GTP');
indexGTP = stateindexIQM(modelfrontend,'G_GTP');
indexGa = stateindexIQM(modelfrontend,'Ga_GTP');

PDE_singcell1 = cell(1,stochsimfrontend.runs);
PDE_singcell2 = cell(1,stochsimfrontend.runs);
PDE_doubcell = cell(1,stochsimfrontend.runs);
timefrontstochALL = cell(1,stochsimfrontend.runs);

G_cell = cell(14, stochsimfrontend.runs);

for k=1:stochsimfrontend.runs
    PDE_sing1 = stochsimfrontend.speciesdata{k}(:,index1);
    PDE_sing2 = stochsimfrontend.speciesdata{k}(:,index3);
    PDE_doub = stochsimfrontend.speciesdata{k}(:,index2);
    timestochk = stochsimfrontend.time{k};
    PDE_singcell1{k} = PDE_sing1;
    PDE_singcell2{k} = PDE_sing2;
    PDE_doubcell{k} = PDE_doub;
    timefrontstochALL{k} = timestochk;
    G_cell{1,k} = stochsimfrontend.speciesdata{k}(:,indexGt0);
    G_cell{2,k} = stochsimfrontend.speciesdata{k}(:,indexGt1);
    G_cell{3,k} = stochsimfrontend.speciesdata{k}(:,indexGt2);
    G_cell{4,k} = stochsimfrontend.speciesdata{k}(:,indexGt3);
    G_cell{5,k} = stochsimfrontend.speciesdata{k}(:,indexG0);
    G_cell{6,k} = stochsimfrontend.speciesdata{k}(:,indexG1);
    G_cell{7,k} = stochsimfrontend.speciesdata{k}(:,indexG2);
    G_cell{8,k} = stochsimfrontend.speciesdata{k}(:,indexG3);
    G_cell{9,k} = stochsimfrontend.speciesdata{k}(:,indexG0GTP);
    G_cell{10,k} = stochsimfrontend.speciesdata{k}(:,indexG1GTP);
    G_cell{11,k} = stochsimfrontend.speciesdata{k}(:,indexG2GTP);
    G_cell{12,k} = stochsimfrontend.speciesdata{k}(:,indexG3GTP);
    G_cell{13,k} = stochsimfrontend.speciesdata{k}(:,indexGTP);
    G_cell{14,k} = stochsimfrontend.speciesdata{k}(:,indexGa);
end
clearvars stochsimfrontend PDE_sing1 PDE_sing2 PDE_doub timestochk modelfrontend;

stochsimfrontend_prec = IQMstochsim2(modelfrontend_prec,V,timeend,runs,Nsample); 

index1 = stateindexIQM(modelfrontend_prec,'PDE_a_Ga_GTP');
index2 = stateindexIQM(modelfrontend_prec,'Ga_GTP_a_PDE_a_Ga_GTP');
index3 = stateindexIQM(modelfrontend_prec, 'Ga_GTP_PDE_a_Ga_GTP');
indexGt0 = stateindexIQM(modelfrontend_prec,'R0_Gt');
indexGt1 = stateindexIQM(modelfrontend_prec,'R1_Gt');
indexGt2 = stateindexIQM(modelfrontend_prec,'R2_Gt');
indexGt3 = stateindexIQM(modelfrontend_prec,'R3_Gt');
indexG0 = stateindexIQM(modelfrontend_prec,'R0_G');
indexG1 = stateindexIQM(modelfrontend_prec,'R1_G');
indexG2 = stateindexIQM(modelfrontend_prec,'R2_G');
indexG3 = stateindexIQM(modelfrontend_prec,'R3_G');
indexG0GTP = stateindexIQM(modelfrontend_prec,'R0_G_GTP');
indexG1GTP = stateindexIQM(modelfrontend_prec,'R1_G_GTP');
indexG2GTP = stateindexIQM(modelfrontend_prec,'R2_G_GTP');
indexG3GTP = stateindexIQM(modelfrontend_prec,'R3_G_GTP');
indexGTP = stateindexIQM(modelfrontend_prec,'G_GTP');
indexGa = stateindexIQM(modelfrontend_prec,'Ga_GTP');

PDE_singcell1_prec = cell(1,stochsimfrontend_prec.runs);
PDE_singcell2_prec = cell(1,stochsimfrontend_prec.runs);
PDE_doubcell_prec = cell(1,stochsimfrontend_prec.runs);
timefrontstochALL_prec = cell(1,stochsimfrontend_prec.runs);

G_cell_prec = cell(14, stochsimfrontend_prec.runs);

for k=1:stochsimfrontend_prec.runs
    PDE_sing1 = stochsimfrontend_prec.speciesdata{k}(:,index1);
    PDE_sing2 = stochsimfrontend_prec.speciesdata{k}(:,index3);
    PDE_doub = stochsimfrontend_prec.speciesdata{k}(:,index2);
    timestochk = stochsimfrontend_prec.time{k};
    PDE_singcell1_prec{k} = PDE_sing1;
    PDE_singcell2_prec{k} = PDE_sing2;
    PDE_doubcell_prec{k} = PDE_doub;
    timefrontstochALL_prec{k} = timestochk;
    G_cell_prec{1,k} = stochsimfrontend_prec.speciesdata{k}(:,indexGt0);
    G_cell_prec{2,k} = stochsimfrontend_prec.speciesdata{k}(:,indexGt1);
    G_cell_prec{3,k} = stochsimfrontend_prec.speciesdata{k}(:,indexGt2);
    G_cell_prec{4,k} = stochsimfrontend_prec.speciesdata{k}(:,indexGt3);
    G_cell_prec{5,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG0);
    G_cell_prec{6,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG1);
    G_cell_prec{7,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG2);
    G_cell_prec{8,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG3);
    G_cell_prec{9,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG0GTP);
    G_cell_prec{10,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG1GTP);
    G_cell_prec{11,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG2GTP);
    G_cell_prec{12,k} = stochsimfrontend_prec.speciesdata{k}(:,indexG3GTP);
    G_cell_prec{13,k} = stochsimfrontend_prec.speciesdata{k}(:,indexGTP);
    G_cell_prec{14,k} = stochsimfrontend_prec.speciesdata{k}(:,indexGa);
end
clearvars stochsimfrontend_prec PDE_sing1 PDE_sing2 PDE_doub timestochk modelfrontend_prec;