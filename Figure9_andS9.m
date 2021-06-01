%% Import data
% 10 data columns
clear all; clc; close all;
fid=fopen('Exp_data/Cell5 combined.axgt');
s=textscan(fid,'%f %f %f %f %f %f %f %f %f %f','headerlines',1);
fclose(fid);

%% plot results
time = s{1,1};
figure(1);
hold on;
plot(time, s{1,2}, 'k');
plot(time, s{1,3}, 'g');
plot(time, s{1,4}, 'r');
plot(time, s{1,5}, 'b');
xlabel('time/s');
ylabel('Green LED/V');

time = s{1,6};

figure(2);
hold on;
plot(time, s{1,7}, 'k');
plot(time, s{1,8}, 'g');
plot(time, s{1,9}, 'r');
plot(time, s{1,10}, 'b');
xlabel('time/s');
ylabel('Photovoltage');

%% Convert to photocurrent
for i=7:10
    [ReCurrent(:,i), ImCurrent(:,i)] = conversionscript(s{1,1}, s{1,i});
end

figure(3);
hold on;
plot(time(1:length(ReCurrent(:,7))), ReCurrent(:,7), 'k');
plot(time(1:length(ReCurrent(:,8))), ReCurrent(:,8), 'g');
plot(time(1:length(ReCurrent(:,9))), ReCurrent(:,9), 'r');
plot(time(1:length(ReCurrent(:,10))), ReCurrent(:,10), 'b');
xlabel('time/s');
ylabel('photocurrent');

%% filter photocurrent
for i=1:4
    current_filt(:,i) = lowpass(ReCurrent(:,i+6),40,5000);
end

%% sliding window average
photocurrent_nice(:,1) = movmean(current_filt(:,1), 100);
photocurrent_nice(:,2) = movmean(current_filt(:,2), 100);
photocurrent_nice(:,3) = movmean(current_filt(:,3), 100);
photocurrent_nice(:,4) = movmean(current_filt(:,4), 100);

%%
figure(4); clf;
hold on;
plot(time(1:length(ReCurrent(:,7))), photocurrent_nice(:,1), 'k');
plot(time(1:length(ReCurrent(:,7))), photocurrent_nice(:,2), 'g');
plot(time(1:length(ReCurrent(:,7))), photocurrent_nice(:,3), 'r');
plot(time(1:length(ReCurrent(:,7))), photocurrent_nice(:,4), 'b');
xlabel('time/s');
ylabel('photocurrent');

%% filter photovoltage
for i=1:4
    voltage_filt(:,i) = lowpass(s{1,i+6},40,5000);
end

%% sliding window average
voltage_nice(:,1) = movmean(voltage_filt(:,1), 100);
voltage_nice(:,2) = movmean(voltage_filt(:,2), 100);
voltage_nice(:,3) = movmean(voltage_filt(:,3), 100);
voltage_nice(:,4) = movmean(voltage_filt(:,4), 100);

%%
figure(5); clf;
hold on;
plot(time, voltage_nice(:,1), 'k');
plot(time, voltage_nice(:,2), 'g');
plot(time, voltage_nice(:,3), 'r');
plot(time, voltage_nice(:,4), 'b');
xlabel('time/s');
ylabel('photovoltage');


%% PUBLICATION FIGURE 9A
figure(6); clf;
hold on;
pbaspect([5 1 1]);
plot(time, voltage_nice(:,1), 'k', 'Linewidth', 1.5);
plot(time, voltage_nice(:,2), 'g', 'Linewidth', 1.5);
plot(time, voltage_nice(:,3), 'r', 'Linewidth', 1.5);
plot(time, voltage_nice(:,4), 'b', 'Linewidth', 1.5);
xlim([4 38]);
ylim([-1.1 0.2])
set(gca, 'Fontsize', 15);
legend('BG 0 ph/{\mu}m^2s', 'BG 8.5 ph/{\mu}m^2s', 'BG 252 ph/{\mu}m^2s', 'BG 695 ph/{\mu}m^2s', 'Location', 'southwest');
xlabel('time/s');
ylabel('scaled {\Delta}U');
text(0.025,0.95,'A','Units','normalized','FontSize',20)

%% simulation
% Recreate light adaptation experiments by Sabrina and Lorenzo

model = IQMmodel('Models/DM.txtbc');
experiment = IQMexperiment('Sim_experiments/flash_background.exp');
modelexp = IQMmergemodexp(model, experiment);

scalingfactor = 0.4;
% old: 0.4
scaling2 = 3;
scalingfactor_bg = 0.4;

flashMags = scalingfactor * [9.5 32.1 96.4 299 5104]; % photons per mum^2
flashDurs = [1e-3 1e-3 3e-3 5e-3 60e-3];
flashBGs = scalingfactor_bg * [0 scaling2*8.5 252 695]; % photons per mum^2 per s

%% adapt model to the different backgrounds
modeladapt = cell(4,1);

% simulate backgrounds
time = [0:0.01:100];
% zero is already dark adapted
modeladapt{1} = modelexp;
modeladapt{1} = IQMparameters(modeladapt{1}, 'flashBG', 0);
modeladapt{1} = IQMparameters(modeladapt{1}, 'flashMag', 0);
% 2-4: adapt to background
for i=2:length(flashBGs)
    modeladapt{i} = IQMparameters(modelexp, 'flashMag', 0);
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashBG', flashBGs(i));
    simulation_adapt(i) = IQMPsimulate(modeladapt{i}, time);
end
% important: for new beelen less E use 76 states
% change starting values in model
for i=2:4
    model_temp = struct(modeladapt{i});
    for j=1:76
        model_temp.states(j).initialCondition = simulation_adapt(i).statevalues(10001,j);
    end
    modeladapt{i} = IQMmodel(model_temp);
end

%% simulate flashes with adapted models

% simulate first flash
time = [0:0.001:3];
for i=1:length(flashBGs)
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashMag', flashMags(1));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashDur', flashDurs(1));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashBG', flashBGs(i));
    simulation1(i) = IQMPsimulate(modeladapt{i}, time);
end

%% simulate second flash
for i=1:length(flashBGs)
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashMag', flashMags(2));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashDur', flashDurs(2));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashBG', flashBGs(i));
    simulation2(i) = IQMPsimulate(modeladapt{i}, time);
end

% simulate third flash
for i=1:length(flashBGs)
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashMag', flashMags(3));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashDur', flashDurs(3));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashBG', flashBGs(i));
    simulation3(i) = IQMPsimulate(modeladapt{i}, time);
end

%% simulate fourth flash
for i=1:length(flashBGs)
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashMag', flashMags(4));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashDur', flashDurs(4));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashBG', flashBGs(i));
    simulation4(i) = IQMPsimulate(modeladapt{i}, time);
end

% simulate fifth flash
for i=1:length(flashBGs)
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashMag', flashMags(5));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashDur', flashDurs(5));
    modeladapt{i} = IQMparameters(modeladapt{i}, 'flashBG', flashBGs(i));
    simulation5(i) = IQMPsimulate(modeladapt{i}, time);
end

%% normalize all plots: first subtract background J
for i = 2:length(flashBGs)
    backgroundJ{i} = simulation_adapt(i).variablevalues(10001, variableindexIQM(modelexp, 'deltaJ'));
end
backgroundJ{1} = 0;

% determine simulation response max to EACH SINGLE FLASH
for i = 1:length(flashBGs)
    maxresponse{1,i} = max(simulation1(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i});
    maxresponse{2,i} = max(simulation2(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i});
    maxresponse{3,i} = max(simulation3(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i});
    maxresponse{4,i} = max(simulation4(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i});
    maxresponse{5,i} = max(simulation5(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i});
end
%% determine measurement response max to EACH SINGLE FLASH
for i = 1:length(flashBGs)
    measmax{1,i} = max(photocurrent_nice(6000:11000,i));
    measmax{2,i} = max(photocurrent_nice(26000:31000,i));
    measmax{3,i} = max(photocurrent_nice(56000:61000,i));
    measmax{4,i} = max(photocurrent_nice(96000:102500,i));
    measmax{5,i} = max(photocurrent_nice(144500:154500,i));
end

%%

colour = ['k', 'g', 'r', 'b'];
for i=1:length(flashBGs)
    deltaJ1{i} = (simulation1(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i})/maxresponse{1,i}*measmax{1,i}/measmax{5,i};
    deltaJ2{i} = (simulation2(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i})/maxresponse{2,i}*measmax{2,i}/measmax{5,i};
    deltaJ3{i} = (simulation3(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i})/maxresponse{3,i}*measmax{3,i}/measmax{5,i};
    deltaJ4{i} = (simulation4(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i})/maxresponse{4,i}*measmax{4,i}/measmax{5,i};
    deltaJ5{i} = (simulation5(i).variablevalues(:, variableindexIQM(modelexp, 'deltaJ'))...
        -backgroundJ{i})/maxresponse{5,i}*measmax{5,i}/measmax{5,i};
end

%% PUBLICATION FIGURE 9B-E
% flash 1 and 2
figure(7);
subplot(2,2,1);
hold on;
text(0.025,0.95,'B','Units','normalized','FontSize',10)
pbaspect([1 1 1]);
plot(s{1,6}(6000:16000), photocurrent_nice(6000:16000,1)/measmax{5,1}, 'k', 'linewidth', 1.5);
plot(s{1,6}(6000:16000), photocurrent_nice(6000:16000,2)/measmax{5,2}, 'g', 'linewidth', 1.5);
plot(s{1,6}(6000:16000), photocurrent_nice(6000:16000,3)/measmax{5,3}, 'r', 'linewidth', 1.5);
plot(s{1,6}(6000:16000), photocurrent_nice(6000:16000,4)/measmax{5,4}, 'b', 'linewidth', 1.5);
ylabel('Scaled {\Delta}J');
xlabel('time/s');
ylim([-0.02 0.12]);
xlim([4.7 6.7]);
set(gca,'Fontsize',12);
title('Flash 1, exp.');
legend('BG 0 ph/{\mu}m^2', 'BG 8.5 ph/{\mu}m^2', 'BG 252 ph/{\mu}m^2', 'BG 695 ph/{\mu}m^2', 'Location', 'northeast');
subplot(2,2,2);
hold on;
text(0.025,0.95,'C','Units','normalized','FontSize',10)
pbaspect([1 1 1]);
plot(time(4700:6700), deltaJall(1,4700:6700), 'k', 'linewidth', 1.5);
plot(time(4700:6700), deltaJall(2,4700:6700), 'g', 'linewidth', 1.5);
plot(time(4700:6700), deltaJall(3,4700:6700), 'r', 'linewidth', 1.5);
plot(time(4700:6700), deltaJall(4,4700:6700), 'b', 'linewidth', 1.5);
ylim([-0.02 0.12]);
xlim([4.7 6.7]);
ylabel('Scaled {\Delta}J');
xlabel('time/s');
set(gca,'Fontsize',12);
title('Flash 1, sim.');

subplot(2,2,3);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'D','Units','normalized','FontSize',10)
plot(s{1,6}(26000:36000), photocurrent_nice(26000:36000,1)/measmax{5,1}, 'k', 'linewidth', 1.5);
plot(s{1,6}(26000:36000), photocurrent_nice(26000:36000,2)/measmax{5,2}, 'g', 'linewidth', 1.5);
plot(s{1,6}(26000:36000), photocurrent_nice(26000:36000,3)/measmax{5,3}, 'r', 'linewidth', 1.5);
plot(s{1,6}(26000:36000), photocurrent_nice(26000:36000,4)/measmax{5,4}, 'b', 'linewidth', 1.5);
ylabel('Scaled {\Delta}J');
xlabel('time/s');
%ylim([-0.02 0.3]);
xlim([8.7 10.7]);
set(gca,'Fontsize',12);
title('Flash 2, exp.');
subplot(2,2,4);
hold on;
text(0.025,0.95,'E','Units','normalized','FontSize',10)
pbaspect([1 1 1]);
plot(time(8700:10700), deltaJall(1,8700:10700), 'k', 'linewidth', 1.5);
plot(time(8700:10700), deltaJall(2,8700:10700), 'g', 'linewidth', 1.5);
plot(time(8700:10700), deltaJall(3,8700:10700), 'r', 'linewidth', 1.5);
plot(time(8700:10700), deltaJall(4,8700:10700), 'b', 'linewidth', 1.5);
%ylim([-0.02 0.3]);
xlim([8.7 10.7]);
ylabel('Scaled {\Delta}J');
xlabel('time/s');
set(gca,'Fontsize',12);
title('Flash 2, sim.');

%% Photovoltage: determine measurement response max to EACH SINGLE FLASH
for i = 1:length(flashBGs)
    measmax_V{1,i} = max(-voltage_nice(6000:11000,i));
    measmax_V{2,i} = max(-voltage_nice(26000:31000,i));
    measmax_V{3,i} = max(-voltage_nice(56000:61000,i));
    measmax_V{4,i} = max(-voltage_nice(96000:102500,i));
    measmax_V{5,i} = max(-voltage_nice(144500:154500,i));
end
%% determine half-maximal width of curve
flashstarts_exp = [6000 26000 56000 96000 144500];
flashstarts_sim = [4900 8900 14900 22900 32900];

for i = 1:4 % all backgrounds
    for j = 1:5 % all flashes
        halfmax_V(j,i) = measmax_V{j,i}/(2*measmax_V{5,i});
        % both simulation and experiment are normalized to this value
        t_exp = flashstarts_exp(j);
        while voltage_nice(t_exp,i)/measmax_V{5,i} > -halfmax_V(j,i)
            t_exp = t_exp+1;
        end
        tend_exp = t_exp+500;
        while voltage_nice(tend_exp,i)/measmax_V{5,i} < -halfmax_V(j,i)
            tend_exp = tend_exp+1;
        end
        % maybe use fit for experimental data because they are so noisy
        halfmax(j,i) = measmax{j,i}/(2*measmax{5,i});
        t_sim = flashstarts_sim(j);
        while deltaJall(i,t_sim) < halfmax(j,i)
            t_sim = t_sim+1;
        end
        tend_sim = t_sim+1;
        while deltaJall(i,tend_sim) > halfmax(j,i)
            tend_sim = tend_sim+1;
        end
        % conversion to time, not data points
        halfmaxtime_exp(j,i) = (tend_exp - t_exp) * 2e-04;
        halfmaxtime_sim(j,i) = (tend_sim - t_sim) * 0.001;
    end
end

%% calculate reduction in halfmax time with respect to dark background time
for i = 1:4 % all backgrounds
    for j = 1:5 % all flashes
        redux_exp(j,i) = halfmaxtime_exp(j,i)/halfmaxtime_exp(j,1);
        redux_sim(j,i) = halfmaxtime_sim(j,i)/halfmaxtime_sim(j,1);
    end
end

%% PUBLICATION FIGURE S9
figure(8);
subplot(1,2,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15)
plot(log(scalingfactor*flashMags), redux_exp(:,1), 'k-o');
plot(log(scalingfactor*flashMags), redux_exp(:,2), 'g-o');
plot(log(scalingfactor*flashMags(:)), redux_exp(:,3), 'r-o');
plot(log(scalingfactor*flashMags(:)), redux_exp(:,4), 'b-o');
ylim([0.4 1.2]);
title('Experiment');
xlabel('log(L) of flash');
ylabel('T_{half} reduction');
set(gca, 'FontSize', 12);

subplot(1,2,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15)
plot(log(scalingfactor*flashMags), redux_sim(:,1), 'k-o');
plot(log(scalingfactor*flashMags), redux_sim(:,2), 'g-o');
plot(log(scalingfactor*flashMags), redux_sim(:,3), 'r-o');
plot(log(scalingfactor*flashMags), redux_sim(:,4), 'b-o');
ylim([0.4 1.2]);
title('Simulation');
xlabel('log(L) of flash');
ylabel('T_{half} reduction');
legend('BG 0 ph/{\mu}m^2s', 'BG 8.5 ph/{\mu}m^2s', 'BG 252 ph/{\mu}m^2s', ...
    'BG 695 ph/{\mu}m^2s', 'Location', 'southeast');
set(gca, 'FontSize', 12);
