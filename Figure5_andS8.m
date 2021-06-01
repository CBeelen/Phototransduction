%% compare stoch sim and electrophys.
clear all; clc;

% load data from stochastic simulation (run on cluster)
load('Simulation_data/Fig5data.mat');
% load electrophysiological data (average SPR from Fig4)
load('Exp_data/avSPR_exp.mat');

%% plot results from stochastic simulations
figure(1); clf;
hold on;
for i=1:num_sim
    plot(time_ALL{1,i}, PDE_doub_ALL{1,i})
end
xlim([0 3]);
xlabel('time');
ylabel('PDE **');
legend('Double activated PDE');

figure(2);
hold on;
for i=1:num_sim
    plot(time_ALL{1,i}, PDE_doub_ALL{1,i}+0.025*PDE_sing_ALL{1,i})
end
xlim([0 3]);
xlabel('time');
ylabel('Effector');
legend('Double activated PDE plus 0.025 single activated');

%% resample PDE
timeend = 3;
sampling = 0.001;
lengthvec = timeend/sampling+1;
PDE_doub_new = zeros(num_sim,lengthvec);
PDE_sing_new = zeros(num_sim,lengthvec);
for i=1:num_sim
    [time_new(:), PDE_doub_new(i,:)] = ...
    resamplePDE(timeend, sampling, time_ALL{1,i}, PDE_doub_ALL{1,i});
    [time_new(:), PDE_sing_new(i,:)] = ...
        resamplePDE(timeend, sampling, time_ALL{1,i}, PDE_sing_ALL{1,i});
end

%% average PDE_doub sing and effector
av_PDE_doub = zeros(1,lengthvec);
for i=1:num_sim
    av_PDE_doub = av_PDE_doub + PDE_doub_new(i,:)/num_sim;
end

av_PDE_sing = zeros(1,lengthvec);
for i=1:num_sim
    av_PDE_sing = av_PDE_sing + PDE_sing_new(i,:)/num_sim;
end

av_effector = av_PDE_doub + 0.025*av_PDE_sing;

%% plot those. PUBLICATION PLOT S8
figure(3); clf;
subplot(1,3,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15);
for i=1:num_sim
    p1 = plot(time_ALL{1,i}, PDE_doub_ALL{1,i}, 'k');
end
p2 = plot(time_new, av_PDE_doub, 'r', 'LineWidth', 1.5);
xlim([0 2]);
ylim([0 60]);
xlabel('time/s');
ylabel('#PDE **');
title('Double activated PDE');
legend([p1 p2], 'Stochastic', 'Average');
set(gca, 'FontSize', 12);

subplot(1,3,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15);
for i=1:num_sim
    p1 = plot(time_ALL{1,i}, PDE_sing_ALL{1,i}, 'k');
end
p2 = plot(time_new, av_PDE_sing, 'r', 'LineWidth', 1.5);
xlim([0 2]);
ylim([0 800]);
xlabel('time/s');
ylabel('#PDE *');
legend([p1 p2], 'Stochastic', 'Average');
title('Single activated PDE');
set(gca, 'FontSize', 12);

subplot(1,3,3);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'C','Units','normalized','FontSize',15);
for i=1:num_sim
    p1 = plot(time_ALL{1,i}, PDE_sing_ALL{1,i} + 2*PDE_doub_ALL{1,i}, 'k');
end
p2 = plot(time_new, av_PDE_sing +2* av_PDE_doub, 'r', 'LineWidth', 1.5);
xlim([0 2]);
ylim([0 900]);
xlabel('time/s');
ylabel('#G *');
legend([p1 p2], 'Stochastic', 'Average');
title('Activated Transducin');
set(gca, 'FontSize', 12);

%% run backend simulation for photocurrent (this takes a few minutes)
modelbackend = IQMmodel('Models/HSDM_backend.txtbc');
deltaJbackend = zeros(301,num_sim);
time = 0:0.01:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL{1,i}),...
    sprintf('%g,',PDE_doub_ALL{1,i}+0.025*PDE_sing_ALL{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% plot results
figure(4);
hold on;
for i=1:num_sim
    plot(time, deltaJbackend(:,i),'r','LineWidth', 1.5);
end
xlim([0 1]);
xlabel('time','Fontsize', 20);
ylabel('\Delta J','Fontsize', 20);
set(gca,'Fontsize',20);


%% calculate average
average = deltaJbackend(:,1)/num_sim;
for i=2:num_sim
    average = average + deltaJbackend(:,i)/num_sim;
end
figure(5);
hold on;
plot(time, average, 'k', 'LineWidth', 1.5);
set(gca,'Fontsize',20);
xlim([0 1]);
xlabel('time/s');
ylabel('\Delta J/pA');

%% calculate coefficient of variability
[CVarea, CVamplitude] = CV(0.01, deltaJbackend(1:150,:), 0,0)

%% compare with deterministic result
model_det = IQMmodel('Models/HSDM_full.txtbc');
sim_det = IQMsimulate(model_det,time);
deltaJ_det = sim_det.variablevalues(:,variableindexIQM(model_det,'deltaJ'));
%% compare stoch and det
figure(6);
hold on;
plot(time, average, 'k', 'LineWidth', 1.5);
plot(time, deltaJ_det, 'r', 'LineWidth', 1.5);
set(gca,'Fontsize',20);
xlim([0 1]);
xlabel('time/s');
ylabel('\Delta J/pA');

%% PUBLICATION FIGURE 5
figure(7);clf;
subplot(1,2,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15);
for i=1:num_sim
    p1 = plot(time, deltaJbackend(:,i),'k','LineWidth', 1);
end
p2 = plot(time, average, 'r','LineWidth', 1.5);
p3 = plot(time, deltaJ_det,'g--','LineWidth', 1.5);
p4 = plot(time1-0.5326, (avSPR-avSPR(2663))/(max(avSPR(2500:5000)-avSPR(2663)))*max(average),...
    'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5);
legend([p1 p2 p3 p4], 'Stochastic', 'Average', 'Deterministic', 'Electrophysiological');
set(gca,'Fontsize',12);
xlim([0 1.5]);
ylim([-0.1 1.7]);
xlabel('time/s');
ylabel('{\Delta}J/pA');

subplot(1,2,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15);
hold on;
for i=1:num_sim
    p1 = plot(time_ALL{1,i}, PDE_doub_ALL{1,i} + 0.025*PDE_sing_ALL{1,i}, 'k');
end
p2 = plot(time_new, av_effector, 'r', 'LineWidth', 1.5);
xlim([0 2]);
ylim([0 70]);
xlabel('time/s');
ylabel('#E');
set(gca,'Fontsize',12);
legend([p1 p2], 'Stochastic', 'Average');

%% save average from stoch sim for other plots
save('Simulation_data/avstochsim.mat', 'average');