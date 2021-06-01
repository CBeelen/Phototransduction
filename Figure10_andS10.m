%% Comparison of model with and without Ca feedback on Rec
clear all; clc; close all;

% integrator options
options = [];
options.maxstep = 0.001;
options.abstol = 1e-8;
options.reltol = 1e-8;

% Load the models

models = [IQMmodel('Models/DM.txtbc'), ...
    IQMmodel('Models/DM_noFB.txtbc')];

mod_cols = {[0 0 0] [0 0 0]};
mod_lcols = ['k' 'k'];
mod_lstyles = {'-' '-'};
mod_symbs = {'o' '*' 's', 'x'};
mod_names = {'DM', 'DM without feedback'};

% Load the experiments
exp_dir = 'Sim_experiments/';
exp_lsf0_less = IQMexperiment([exp_dir 'light_step_and_flash_0_less.exp']);
exp_lsf1_less = IQMexperiment([exp_dir 'light_step_and_flash_1_less.exp']);
exp_lsf2_less = IQMexperiment([exp_dir 'light_step_and_flash_2_less.exp']);
exp_lsf3_less = IQMexperiment([exp_dir 'light_step_and_flash_3_less.exp']);

models_WT_lsf0_less = [];
models_WT_lsf1_less = [];
models_WT_lsf2_less = [];
models_WT_lsf3_less = [];

for m=1:length(models)
  models_WT_lsf0_less = [models_WT_lsf0_less IQMmergemodexp(models(m), exp_lsf0_less)];
  models_WT_lsf1_less = [models_WT_lsf1_less IQMmergemodexp(models(m), exp_lsf1_less)];
  models_WT_lsf2_less = [models_WT_lsf2_less IQMmergemodexp(models(m), exp_lsf2_less)];
  models_WT_lsf3_less = [models_WT_lsf3_less IQMmergemodexp(models(m), exp_lsf3_less)];
end
          
%% Simulation

timestep = 0.001;
for m=1:length(models)
    y0_less = IQMPsimulate(models_WT_lsf0_less(m), [0:timestep:30]);
    deltaJ0_less{m} = y0_less.variablevalues(:, variableindexIQM(models_WT_lsf0_less(m), 'deltaJ'));
    y1_less = IQMPsimulate(models_WT_lsf1_less(m), [0:timestep:30]);
    deltaJ1_less{m} = y1_less.variablevalues(:, variableindexIQM(models_WT_lsf1_less(m), 'deltaJ'));
    y2_less = IQMPsimulate(models_WT_lsf2_less(m), [0:timestep:30]);
    deltaJ2_less{m} = y2_less.variablevalues(:, variableindexIQM(models_WT_lsf2_less(m), 'deltaJ'));
    y3_less = IQMPsimulate(models_WT_lsf3_less(m), [0:timestep:30]);
    deltaJ3_less{m} = y3_less.variablevalues(:, variableindexIQM(models_WT_lsf3_less(m), 'deltaJ'));
end

%% PUBLICATION FIGURE 10
figure(1); clf;
subplot(1,2,1);
hold on;
pbaspect([2 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15);
xlim([9.5 13]);
ylim([-2 16]);
xlabel('time/s');
ylabel('{\Delta}J/pA');
set(gca,'FontSize',15)
plot_lsf = plot(y0_less.time, deltaJ0_less{1}, '-', 'LineWidth', 2, 'Color', 'k');
plot_lsf1 = plot(y1_less.time, deltaJ1_less{1}, '--', 'LineWidth', 2, 'Color', 'g');
plot_lsf2 = plot(y2_less.time, deltaJ2_less{1}, '-.', 'LineWidth', 2, 'Color', 'r');
plot_lsf3 = plot(y3_less.time, deltaJ3_less{1}, ':', 'LineWidth', 2, 'Color', 'b');
title('Normal model');

subplot(1,2,2);
hold on;
pbaspect([2 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15);
xlim([9.5 13]);
ylim([-2 16]);
xlabel('time/s');
ylabel('{\Delta}J/pA');
set(gca,'FontSize',15)
plot_lsf = plot(y0_less.time, deltaJ0_less{2}, '-', 'LineWidth', 2, 'Color', 'k');
plot_lsf1 = plot(y1_less.time, deltaJ1_less{2}, '--', 'LineWidth', 2, 'Color', 'g');
plot_lsf2 = plot(y2_less.time, deltaJ2_less{2}, '-.', 'LineWidth', 2, 'Color', 'r');
plot_lsf3 = plot(y3_less.time, deltaJ3_less{2}, ':', 'LineWidth', 2, 'Color', 'b');
legend('BG 0 ph/{\mu}m^2s', 'BG 698 ph/{\mu}m^2s', 'BG 1860 ph/{\mu}m^2s', ...
    'BG 4651 ph/{\mu}m^2s', 'Location', 'northeast');
title('Calcium Feedback on Rec missing');

%% more data points for saturation times
% Load the experiments
exp_lsf0 = IQMexperiment([exp_dir 'light_step_and_flash_0.exp']);
exp_lsf1 = IQMexperiment([exp_dir 'light_step_and_flash_1.exp']);
exp_lsf2 = IQMexperiment([exp_dir 'light_step_and_flash_2.exp']);
exp_lsf3 = IQMexperiment([exp_dir 'light_step_and_flash_3.exp']);
exp_lsf4 = IQMexperiment([exp_dir 'light_step_and_flash_4.exp']);
exp_lsf5 = IQMexperiment([exp_dir 'light_step_and_flash_5.exp']);
exp_lsf6 = IQMexperiment([exp_dir 'light_step_and_flash_6.exp']);

models_WT_lsf0 = [];
models_WT_lsf1 = [];
models_WT_lsf2 = [];
models_WT_lsf3 = [];
models_WT_lsf4 = [];
models_WT_lsf5 = [];
models_WT_lsf6 = [];
for m=1:length(models)
  models_WT_lsf0 = [models_WT_lsf0 IQMmergemodexp(models(m), exp_lsf0)];
  models_WT_lsf1 = [models_WT_lsf1 IQMmergemodexp(models(m), exp_lsf1)];
  models_WT_lsf2 = [models_WT_lsf2 IQMmergemodexp(models(m), exp_lsf2)];
  models_WT_lsf3 = [models_WT_lsf3 IQMmergemodexp(models(m), exp_lsf3)];
  models_WT_lsf4 = [models_WT_lsf4 IQMmergemodexp(models(m), exp_lsf4)];
  models_WT_lsf5 = [models_WT_lsf5 IQMmergemodexp(models(m), exp_lsf5)];
  models_WT_lsf6 = [models_WT_lsf6 IQMmergemodexp(models(m), exp_lsf6)];
end

%% Simulation

timestep = 0.001;
for m=1:length(models)
    y0 = IQMPsimulate(models_WT_lsf0(m), [0:timestep:30]);
    deltaJ0{m} = y0.variablevalues(:, variableindexIQM(models_WT_lsf0(m), 'deltaJ'));
    y1 = IQMPsimulate(models_WT_lsf1(m), [0:timestep:30]);
    deltaJ1{m} = y1.variablevalues(:, variableindexIQM(models_WT_lsf1(m), 'deltaJ'));
    y2 = IQMPsimulate(models_WT_lsf2(m), [0:timestep:30]);
    deltaJ2{m} = y2.variablevalues(:, variableindexIQM(models_WT_lsf2(m), 'deltaJ'));
    y3 = IQMPsimulate(models_WT_lsf3(m), [0:timestep:30]);
    deltaJ3{m} = y3.variablevalues(:, variableindexIQM(models_WT_lsf3(m), 'deltaJ'));
    y4 = IQMPsimulate(models_WT_lsf4(m), [0:timestep:30]);
    deltaJ4{m} = y4.variablevalues(:, variableindexIQM(models_WT_lsf4(m), 'deltaJ'));
    y5 = IQMPsimulate(models_WT_lsf5(m), [0:timestep:30]);
    deltaJ5{m} = y5.variablevalues(:, variableindexIQM(models_WT_lsf5(m), 'deltaJ'));
    y6 = IQMPsimulate(models_WT_lsf6(m), [0:timestep:30]);
    deltaJ6{m} = y6.variablevalues(:, variableindexIQM(models_WT_lsf6(m), 'deltaJ'));
end

%% Calculate T_sat

timeinsat = zeros(7, length(models));

for m=1:length(models)
    amp = max(deltaJ0{m});
    vec_timeinsat = deltaJ0{m} > 0.95*amp;
    timeinsat(1, m) = timestep * sum(vec_timeinsat);
    amp = max(deltaJ1{m});
    vec_timeinsat = deltaJ1{m} > 0.95*amp;
    timeinsat(2,m) = timestep * sum(vec_timeinsat);
    amp = max(deltaJ2{m});
    vec_timeinsat = deltaJ2{m} > 0.95*amp;
    timeinsat(3,m) = timestep * sum(vec_timeinsat);
    amp = max(deltaJ3{m});
    vec_timeinsat = deltaJ3{m} > 0.95*amp;
    timeinsat(4,m) = timestep * sum(vec_timeinsat);
    amp = max(deltaJ4{m});
    vec_timeinsat = deltaJ4{m} > 0.95*amp;
    timeinsat(5,m) = timestep * sum(vec_timeinsat);
    amp = max(deltaJ5{m});
    vec_timeinsat = deltaJ5{m} > 0.95*amp;
    timeinsat(6,m) = timestep * sum(vec_timeinsat);
    amp = max(deltaJ6{m});
    vec_timeinsat = deltaJ6{m} > 0.95*amp;
    timeinsat(7,m) = timestep * sum(vec_timeinsat);
end
flashintensities=[0 300 600 2000 100 1000 1500];

%%
figure(2);
hold on;
plot(log(flashintensities/0.43), timeinsat(:,1), '+', 'Color', 'k', 'MarkerSize', 10);
plot(log(flashintensities/0.43), timeinsat(:,2), '*', 'Color', 'r', 'MarkerSize', 10);
legend('Normal model', 'Calcium Feedback on Rec missing', 'Location', 'southwest');
xlabel('ln(L)','FontSize', 15);
ylabel('T_{sat}/s','FontSize', 15);
set(gca,'FontSize',15);

%% comparison: calcium feedback makes no difference withput background
% Load the experimental settings
exp_kolesnikov = IQMexperiment([exp_dir 'kolesnikov2010_young_s1.exp']);
model_WT_kol = IQMmergemodexp(models(1), exp_kolesnikov);
model_new_kol = IQMmergemodexp(models(2), exp_kolesnikov);
flashMag   = [0.731 2.064 6.536 16.942 53.75 190.92 604.58 1990.9];
deltaJ = zeros(1001, length(flashMag));
deltaJ_new = zeros(1001, length(flashMag));
for k=1:length(flashMag), k
    model_WT_kol_k = IQMparameters(model_WT_kol,'flashMag',flashMag(k));
    output_sim = IQMPsimulate(model_WT_kol_k,[0:0.005:5],[],[],[],options);
    deltaJ(:,k) = output_sim.variablevalues(:,variableindexIQM(model_WT_kol_k,'deltaJ'));
    output_sim.variablevalues(1,variableindexIQM(model_WT_kol_k,'mag'))
    model_new_kol_k = IQMparameters(model_new_kol,'flashMag',flashMag(k));
    output_sim_new = IQMPsimulate(model_new_kol_k,[0:0.005:5],[],[],[],options);
    deltaJ_new(:,k) = output_sim_new.variablevalues(:,variableindexIQM(model_new_kol_k,'deltaJ'));
    output_sim_new.variablevalues(1,variableindexIQM(model_new_kol_k,'mag'))
end
%%
figure(3); clf;
hold on;
for k=1:length(flashMag)
    p1 = plot(output_sim.time,deltaJ(:,k),'k-','LineWidth',1);
    p2 = plot(output_sim_new.time,deltaJ_new(:,k),'r--','LineWidth',1);
end
xlabel('time (s)','FontSize',14,'FontWeight','bold');
ylabel('\DeltaJ (pA)','FontSize',14,'FontWeight','bold');
title('Jarvinen Flash Series','FontWeight','bold');
xlim([-0.25 4])
ylim([-1 16])
legend([p1 p2], 'Normal model', 'Calcium Feedback on Rec missing');
set(gca, 'FontSize', 15);

%% PUBLICATION FIGURE S10
figure(4); clf;
subplot(1,2,1);
hold on;
text(0.025,0.95,'A','Units','normalized','FontSize',18)
pbaspect([1 1 1]);
for k=1:length(flashMag)
    p1 = plot(output_sim.time,deltaJ(:,k),'k-','LineWidth',1.5);
    p2 = plot(output_sim_new.time,deltaJ_new(:,k),'r--','LineWidth',1.5);
end
xlabel('time/s');
ylabel('\DeltaJ/pA');
title('Flash Series without background');
xlim([-0.25 5])
ylim([-1 16])
legend([p1 p2], 'Normal model', 'No Calcium Feedback on Rec');
set(gca,'FontSize',15);

subplot(1,2,2);
hold on;
text(0.025,0.95,'B','Units','normalized','FontSize',18)
pbaspect([1 1 1]);
plot(log(flashintensities/0.43), timeinsat(:,1), '+', 'Color', 'k', 'MarkerSize', 10);
plot(log(flashintensities/0.43), timeinsat(:,2), '*', 'Color', 'r', 'MarkerSize', 10);
legend('Normal model', 'No Calcium Feedback on Rec', 'Location', 'southwest');
xlabel('ln(L)','FontSize', 15);
ylabel('T_{sat}/s','FontSize', 15);
set(gca,'FontSize',15);
title('Saturation times with background');


