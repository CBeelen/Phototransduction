%% Comparison of DM and Invergo models
clear all; clc; close all;

exp_dir = 'Sim_experiments/';

% Load the models
model_WT = IQMmodel('Models/Invergo.txt');
model_new = IQMmodel('Models/DM.txtbc');

% Load the experimental settings
exp_kolesnikov = IQMexperiment([exp_dir 'kolesnikov2010_young_s1.exp']);

% Simulator options
options = [];
options.maxstep = 0.001;
options.abstol = 1e-8;
options.reltol = 1e-8;

% Merge experiment and model
model_WT_kol = IQMmergemodexp(model_WT, exp_kolesnikov);
model_new_kol = IQMmergemodexp(model_new, exp_kolesnikov);
flashMag   = [0.731 2.064 6.536 16.942 53.75 190.92 604.58 1990.9];
deltaJ = zeros(1001, length(flashMag));
deltaJ_new = zeros(1001, length(flashMag));

%% Simulate
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
figure(1); clf;
hold on;
for k=1:length(flashMag)
    p1 = plot(output_sim.time,deltaJ(:,k),'k-','LineWidth',1);
    p2 = plot(output_sim_new.time,deltaJ_new(:,k),'r--','LineWidth',1);
end
xlabel('time (s)','FontSize',14,'FontWeight','bold');
ylabel('\DeltaJ (pA)','FontSize',14,'FontWeight','bold');
title(' Jarvinen Flash Series from 1.54 up to 118,000 R*/flash','FontWeight','bold');
xlim([-0.25 5])
ylim([-1 16])
legend([p1 p2], 'Invergo 2014', 'DM');

%%
% Load the experiments
exp_bg = IQMexperiment([exp_dir 'LIGHT_ADAPT_lightbckg.exp']);
exp_bg_flash = IQMexperiment([exp_dir 'LIGHT_ADAPT_lightbckg_flash.exp']);
exp_flash = IQMexperiment([exp_dir 'LIGHT_ADAPT_onlyflash.exp']);
exp_GCAPs_ko = IQMexperiment([exp_dir 'GCAPs_ko.exp']);
exp_GCAPs_ko_DM = IQMexperiment([exp_dir 'GCAPs_ko_forDM.exp']);

models = [model_WT model_new];
mod_cols = {[0 0 0.5] [0.8 0.8 0.8] [0.5 0 0] [0 0.5 0]};
mod_lcols = ['k' 'r' 'b' 'g'];

models_WT_bg = [];
models_WT_bg_flash = [];
models_WT_flash = [];
models_GCAPs_ko_bg = [];
models_GCAPs_ko_bg_flash = [];
models_GCAPs_ko_flash = [];
% GCAPs KO needs a different experiment for each model, because of the
% changes to the GC-GCAP implementation
model_GCAPs_ko = [IQMmergemodexp(models(1), exp_GCAPs_ko) IQMmergemodexp(models(2), exp_GCAPs_ko_DM)];

for m=1:length(models),
  models_WT_bg = [models_WT_bg IQMmergemodexp(models(m), exp_bg)];
  models_WT_bg_flash = [models_WT_bg_flash IQMmergemodexp(models(m),...
							 exp_bg_flash)];
  models_WT_flash = [models_WT_flash IQMmergemodexp(models(m), exp_flash)];
  models_GCAPs_ko_bg = [models_GCAPs_ko_bg IQMmergemodexp(model_GCAPs_ko(m), exp_bg)];
  models_GCAPs_ko_bg_flash = [models_GCAPs_ko_bg_flash IQMmergemodexp(model_GCAPs_ko(m),...
							       exp_bg_flash)];
  models_GCAPs_ko_flash = [models_GCAPs_ko_flash IQMmergemodexp(model_GCAPs_ko(m),...
							       exp_flash)];
end

%%
deltaJ_bg = zeros(110001,length(models));
deltaJ_bg_flash = zeros(110001,length(models));
deltaJ_flash = zeros(110001,length(models));

for m=1:length(models)
    y1 = IQMPsimulate(models_WT_bg(m), [0:0.001:110]);
    deltaJ_bg(:,m) = y1.variablevalues(:, variableindexIQM(models_WT_bg(m), 'deltaJ'));

    y2 = IQMPsimulate(models_WT_bg_flash(m), [0:0.001:110]);
    deltaJ_bg_flash(:,m) = y2.variablevalues(:, variableindexIQM(models_WT_bg_flash(m), 'deltaJ'));

    y3 = IQMPsimulate(models_WT_flash(m), [0:0.001:110]);
    deltaJ_flash(:,m) = y3.variablevalues(:, variableindexIQM(models_WT_flash(m), 'deltaJ'));
end

deltaJ_bg_GCAP = zeros(110001,length(models));
deltaJ_bg_flash_GCAP = zeros(110001,length(models));
deltaJ_flash_GCAP = zeros(110001,length(models));

for m=1:length(models)
   y1 = IQMPsimulate(models_GCAPs_ko_bg(m), [0:0.001:110]);
   deltaJ_bg_GCAP(:,m) = y1.variablevalues(:, variableindexIQM(models_GCAPs_ko_bg(m), 'deltaJ'));

   y2 = IQMPsimulate(models_GCAPs_ko_bg_flash(m), [0:0.001:110]);
   deltaJ_bg_flash_GCAP(:,m) = y2.variablevalues(:, variableindexIQM(models_GCAPs_ko_bg_flash(m), 'deltaJ'));

   y3 = IQMPsimulate(models_GCAPs_ko_flash(m), [0:0.001:110]);
   deltaJ_flash_GCAP(:,m) = y3.variablevalues(:, variableindexIQM(models_GCAPs_ko_flash(m), 'deltaJ'));
end

%%
figure(2);
subplot(1,2,1);
hold on;
pbaspect([1,1,1]);
for m=1:length(models)
    plot_bg = plot(y1.time, deltaJ_bg(:,m), [mod_lcols(m) '-.'], 'LineWidth', 1);
    plot_bg_flash = plot(y2.time, deltaJ_bg_flash(:,m), [mod_lcols(m) '--'], 'LineWidth', 1);
    plot_flash(m) = plot(y3.time, deltaJ_flash(:,m), [mod_lcols(m) '-'], 'LineWidth', 1);
end
xlim([100 102]);
ylim([4 16]);
xlabel('time (s)','FontWeight','bold','FontSize', 11);
ylabel('Photocurrent (pA)','FontWeight','bold','FontSize', 11);
title('WT');
hold off;
%legend([plot_flash(1) plot_flash(2)], 'Invergo 2014', 'Beelen 2020', 'Location', 'Northeast');

subplot(1,2,2);
hold on;
pbaspect([1,1,1]);
for m=1:length(models)
    plot_bg_GCAP(m) = plot(y1.time, deltaJ_bg_GCAP(:,m), [mod_lcols(m) '-.'], 'LineWidth', 1);
    plot_bg_flash_GCAP(m) = plot(y2.time, deltaJ_bg_flash_GCAP(:,m), [mod_lcols(m) '--'], 'LineWidth', 1);
    plot_flash_GCAP(m) = plot(y3.time, deltaJ_flash_GCAP(:,m), [mod_lcols(m) '-'], 'LineWidth', 1);
end
xlim([100 102]);
ylim([4 16]);
xlabel('time (s)','FontWeight','bold','FontSize', 11);
ylabel('Photocurrent (pA)','FontWeight','bold','FontSize', 11);
title('GCAPs-/-');
legend([plot_bg_GCAP(1) plot_bg_flash_GCAP(1) plot_flash_GCAP(1) plot_flash(1) plot_flash(2)], ...
    'Background','Backgound+flash', 'Flash', 'Invergo 2014', 'Beelen 2020', 'Location', 'Southwest');

%% PUBLICATION FIGURE S6
figure(3); clf;
subplot(1,3,1);
hold on;
pbaspect([1,1,1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15);
for k=1:length(flashMag)
    p1 = plot(output_sim.time,deltaJ(:,k),'k-','LineWidth',1);
    p2 = plot(output_sim_new.time,deltaJ_new(:,k),'r-','LineWidth',1);
end
xlabel('time (s)');
ylabel('Photocurrent (pA)');
title(' Jarvinen Flash Series');
xlim([0 4])
ylim([0 17])
legend([p1 p2], 'Invergo 2014', 'Beelen 2020');
set(gca, 'FontSize', 12)

subplot(1,3,2);
hold on;
pbaspect([1,1,1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15);
for m=1:length(models)
    plot_bg = plot(y1.time, deltaJ_bg(:,m), [mod_lcols(m) '-.'], 'LineWidth', 1);
    plot_bg_flash = plot(y2.time, deltaJ_bg_flash(:,m), [mod_lcols(m) '--'], 'LineWidth', 1);
    plot_flash(m) = plot(y3.time, deltaJ_flash(:,m), [mod_lcols(m) '-'], 'LineWidth', 1);
end
xlim([100 102]);
ylim([4 16]);
xlabel('time (s)');
ylabel('Photocurrent (pA)');
title('Light Adaptation: WT');
set(gca, 'FontSize', 12)

subplot(1,3,3);
hold on;
pbaspect([1,1,1]);
text(0.025,0.95,'C','Units','normalized','FontSize',15);
for m=1:length(models)
    plot_bg_GCAP(m) = plot(y1.time, deltaJ_bg_GCAP(:,m), [mod_lcols(m) '-.'], 'LineWidth', 1);
    plot_bg_flash_GCAP(m) = plot(y2.time, deltaJ_bg_flash_GCAP(:,m), [mod_lcols(m) '--'], 'LineWidth', 1);
    plot_flash_GCAP(m) = plot(y3.time, deltaJ_flash_GCAP(:,m), [mod_lcols(m) '-'], 'LineWidth', 1);
end
xlim([100 102]);
ylim([8 16]);
xlabel('time (s)');
ylabel('Photocurrent (pA)');
title('Light Adaptation: GCAPs-/-');
legend([plot_bg_GCAP(1) plot_bg_flash_GCAP(1) plot_flash_GCAP(1)], ...
    'Background','Background+flash', 'Flash', 'Location', 'Southwest');
set(gca, 'FontSize', 12)



