%% compare experimental data and stoch sim: KO experiments
clear all; clc;

% load data:
load('Simulation_data/Fig8data.mat');
% load stoch sim average SPR:
load('Simulation_data/avstochsim.mat');

num_sim = 100;

%% run backend simulation to get photocurrent: Arrestin K/O
modelbackend = IQMmodel('Models/HSDM_backend.txtbc');
deltaJbackend_Arr = zeros(301,num_sim);
time = 0:0.01:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL_Arr{1,i}),...
    sprintf('%g,',PDE_doub_ALL_Arr{1,i}+0.025*PDE_sing_ALL_Arr{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend_Arr(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% run backend simulation to get photocurrent: RK K/O

modelbackend = IQMmodel('Models/HSDM_backend.txtbc');
deltaJbackend_RK = zeros(301,num_sim);
time = 0:0.01:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL_RK{1,i}),...
    sprintf('%g,',PDE_doub_ALL_RK{1,i}+0.025*PDE_sing_ALL_RK{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend_RK(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% run backend simulation to get photocurrent: CSM

modelbackend = IQMmodel('Models/HSDM_backend.txtbc');
deltaJbackend_CSM = zeros(301,num_sim);
time = 0:0.01:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL_CSM{1,i}),...
    sprintf('%g,',PDE_doub_ALL_CSM{1,i}+0.025*PDE_sing_ALL_CSM{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend_CSM(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% run backend simulation to get photocurrent: GCAPs K/O

modelbackend = IQMmodel('Models/HSDM_backend_GCAPsKO.txtbc');
deltaJbackend_GCAP = zeros(301,num_sim);
time = 0:0.01:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_ALL_GCAP{1,i}),...
    sprintf('%g,',PDE_doub_ALL_GCAP{1,i}+0.025*PDE_sing_ALL_GCAP{1,i}));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJbackend_GCAP(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end


%% calculate averages
average_Arr = deltaJbackend_Arr(:,1)/num_sim;
for i=2:num_sim
    average_Arr = average_Arr + deltaJbackend_Arr(:,i)/num_sim;
end

average_RK = deltaJbackend_RK(:,1)/num_sim;
for i=2:num_sim
    average_RK = average_RK + deltaJbackend_RK(:,i)/num_sim;
end

average_CSM = deltaJbackend_CSM(:,1)/num_sim;
for i=2:num_sim
    average_CSM = average_CSM + deltaJbackend_CSM(:,i)/num_sim;
end

average_GCAP = deltaJbackend_GCAP(:,1)/num_sim;
for i=2:num_sim
    average_GCAP = average_GCAP + deltaJbackend_GCAP(:,i)/num_sim;
end

%% load Grabit data: grabbed from references
load('Data_references/Arr.mat');
load('Data_references/WT_Arr.mat');
load('Data_references/CSM.mat');
load('Data_references/WT_CSM.mat');
load('Data_references/RK.mat');
load('Data_references/WT_RK.mat');
load('Data_references/GCAPs.mat');
load('Data_references/WT_GCAPs.mat');

%% PUBLICATION FIGURE 3
figure(4); clf;
subplot(1,2,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15);
% mutants
plot(Arr(:,1), Arr(:,2)/max(WT_Arr(:,2)), 'r', 'LineWidth', 1);
plot(RK(:,1), RK(:,2)/max(WT_RK(:,2)), 'b', 'LineWidth', 1);
plot(CSM(:,1), CSM(:,2)/max(WT_CSM(:,2)), 'g', 'LineWidth', 1);
plot(GCAPs(:,1), GCAPs(:,2)/max(WT_GCAPs(:,2)), 'Color', ...
    [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
% wild types
plot(WT_Arr(:,1), WT_Arr(:,2)/max(WT_Arr(:,2)), 'r.', 'LineWidth', 0.5);
plot(WT_RK(:,1), WT_RK(:,2)/max(WT_RK(:,2)), 'b.', 'LineWidth', 0.5);
plot(WT_CSM(:,1), WT_CSM(:,2)/max(WT_CSM(:,2)), 'g.', 'LineWidth', 0.5);
plot(WT_GCAPs(:,1), WT_GCAPs(:,2)/max(WT_GCAPs(:,2)), '.', 'Color', ...
    [0.9290, 0.6940, 0.1250], 'LineWidth', 0.5);
xlim([0 1]);
ylim([0 4.5]);
xlabel('time (s)');
ylabel('Scaled Photocurrent');
set(gca,'Fontsize',12);
subplot(1,2,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15);
plot(time, average/max(average), 'k.', 'LineWidth', 1);
plot(time, average_Arr/max(average), 'r', 'LineWidth', 1);
plot(time, average_RK/max(average), 'b', 'LineWidth', 1);
plot(time, average_CSM/max(average), 'g', 'LineWidth', 1);
plot(time, average_GCAP/max(average), 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
xlim([0 1]);
xlabel('time (s)');
ylabel('Scaled Photocurrent');
L = legend('WT', 'Arr -/-', 'Rk -/-', 'CSM', 'GCAPs -/-');
L.FontSize = 10;
L.ItemTokenSize = [15,9];
set(gca,'Fontsize',12);
ylim([0 4]);
