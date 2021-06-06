%% compare precoupled and non-precoupled SPRs
clear all; close all;
% load data:
load('Simulation_data/Fig7data.mat');
num_sim = 1000;

%% backend simulation to get photocurrent (this will take a while!)
modelbackend = IQMmodel('Models/HSDM_backend.txtbc');
deltaJbackend = zeros(301,num_sim);
time = 0:0.01:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL{1,i}),...
        sprintf('%g,',PDE_doub_ALL{1,i}+0.025*PDE_sing1_ALL{1,i}+...
        0.025*PDE_sing2_ALL{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% delta J for precoupled case (this will take a while!)
deltaJbackend_prec = zeros(301,num_sim);
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL_prec{1,i}),...
        sprintf('%g,',PDE_doub_ALL_prec{1,i}+0.025*PDE_sing1_ALL_prec{1,i}+...
        0.025*PDE_sing2_ALL_prec{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend_prec(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% calculate average
average = deltaJbackend(:,1)/num_sim;
for i=2:num_sim
    average = average + deltaJbackend(:,i)/num_sim;
end

average_prec = deltaJbackend_prec(:,1)/num_sim;
for i=2:num_sim
    average_prec = average_prec + deltaJbackend_prec(:,i)/num_sim;
end


%% plot
figure(1);clf;
hold on;
plot(time, average, 'k', 'LineWidth', 1);
plot(time, average_prec, 'r--', 'LineWidth', 1);
set(gca,'Fontsize',15);
xlim([0 1]);
ylim([-0.05 0.7]);
xlabel('time/s');
ylabel('{\Delta}J/pA');
legend('Rh not precoupled', 'Rh precoupled', 'Location', 'Northeast');
hold off;

% create smaller axes in top right, and plot on it
axes('Position',[0.57 0.45 .3 .3]);
box on;
hold on;
plot(time, average, 'k', 'LineWidth', 1);
plot(time, average_prec, 'r--', 'LineWidth', 1);
xlim([0.18,0.28]);
ylim([0.61 0.65]);
hold off;

%% calculate amplitudes
amps = zeros(1, num_sim);
amps_prec = zeros(1, num_sim);
for i = 1:num_sim
    amps(i) = max(deltaJbackend(:,i));
    amps_prec(i) = max(deltaJbackend_prec(:,i));
end

%% calculate coefficients of variability
[CVarea, CVamplitude] = CV(0.01, deltaJbackend, 0,0)
[CVarea_prec, CVamplitude_prec] = CV(0.01, deltaJbackend_prec, 0,0)
stddev = sqrt(var(amps));
stddev_prec = sqrt(var(amps_prec));

%% time to peak
maxindex = zeros(1,num_sim);
maxindex_prec = zeros(1,num_sim);

maxdJ = zeros(1,num_sim);
maxdJ_prec = zeros(1,num_sim);

for i = 1:num_sim
    [maxdJ(i), maxindex(i)] = max(deltaJbackend(:,i));
    [maxdJ_prec(i), maxindex_prec(i)] = max(deltaJbackend_prec(:,i));
end
% convert to times
maxindex = maxindex * 0.01;
maxindex_prec = maxindex_prec * 0.01;

% compute average and stdev
av_ttp = mean(maxindex);
av_ttp_prec = mean(maxindex_prec);
stddev_ttp = sqrt(var(maxindex));
stddev_ttp_prec = sqrt(var(maxindex_prec));
av_amp = mean(maxdJ);
av_amp_prec = mean(maxdJ_prec);
stddev_amp = sqrt(var(maxdJ));
stddev_amp_prec = sqrt(var(maxdJ_prec));

%% calculate areas
timestep = 0.01;
areas = zeros(1, num_sim);
areas_prec = zeros(1, num_sim);
for i= 1:num_sim
    areas(i) = sum(deltaJbackend(:,i)*timestep);
    areas_prec(i) = sum(deltaJbackend_prec(:,i)*timestep);
end

%% PUBLICATION FIGURE 7
edges = 0:0.04:2;
figure(2); clf;
subplot(2,2,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',12);
text(0.5,0.8,'CV = 0.90','Units','normalized','FontSize',12);
histogram(areas, edges);
xlabel('Area (pA*s)');
ylabel('Counts');
title('No pre-assembly');
ylim([0 145]);
subplot(2,2,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',12);
text(0.5,0.8,'CV = 0.85','Units','normalized','FontSize',12);
histogram(areas_prec, edges);
xlabel('Area (pA*s)');
ylabel('Counts');
title('Pre-assembly');
ylim([0 145]);

edges = 0:0.04:1.5;
subplot(2,2,3);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'C','Units','normalized','FontSize',12);
text(0.1,0.8,'CV = 0.36','Units','normalized','FontSize',12);
histogram(amps, edges);
xlabel('Amplitude (pA)');
ylabel('Counts');
title('No pre-assembly');
ylim([0 90]);
subplot(2,2,4);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'D','Units','normalized','FontSize',12);
text(0.1,0.8,'CV = 0.34','Units','normalized','FontSize',12);
histogram(amps_prec, edges);
xlabel('Amplitude (pA)');
ylabel('Counts');
title('Pre-assembly ');
ylim([0 90]);
