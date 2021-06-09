%% load models and exp trace
clc; clear all; close all;
model_small = IQMmodel('Models/HSDM_full.txtbc');
model_test = IQMmodel('Models/HSDM_test.txtbc');

load('Exp_data/avSPR_exp.mat');

%% simulate
time = 0:0.001:3;
sim_small = IQMPsimulate(model_small,time);
deltaJ_small = sim_small.variablevalues(:,variableindexIQM(model_small,'deltaJ'));

sim_test = IQMPsimulate(model_test,time);
deltaJ_test = sim_test.variablevalues(:,variableindexIQM(model_test, 'deltaJ'));

%% PUBLICATION FIGURE S1
figure(3); clf,
hold on;
plot(time, deltaJ_small/max(deltaJ_small), 'k', 'LineWidth', 1.5);
plot(time, deltaJ_test/max(deltaJ_test), 'r--', 'LineWidth', 1.5);
plot(time1-0.533, (avSPR-avSPR(2663))/max(avSPR(2500:5000)-avSPR(2663)), 'c-.', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('Scaled Photocurrent');
xlim([0 1.2]);
ylim([-0.03 1.05]);
legend('HSDM', 'Test', 'Experimental', 'Location', 'Northeast');
set(gca, 'FontSize', 18);

%% TTP
[maxvalue, index] = max(deltaJ_small)
maxPDEdoub_small = max(sim_small.statevalues(:,stateindexIQM(model_small, 'Ga_GTP_a_PDE_a_Ga_GTP')))
[maxvalue_test, index_test] = max(deltaJ_test)
maxPDEdoub_test = max(sim_test.statevalues(:,stateindexIQM(model_test, 'Ga_GTP_a_PDE_a_Ga_GTP')))

