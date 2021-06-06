%% Random initial conditions
clear all; close all;
% load data:
load('Simulation_data/Fig6data.mat');

num_sim = 200;
max_samp = size(time,2);

%% rename
num_prec = G_R0;
num_nonprec = R0;

%% plot PDE
figure(1);clf;
hold on;
for i=1:num_sim
    for j=1:max_samp
        plot(time{i,j}, PDE_doub{i,j}, 'r');
    end
end
xlim([0 3]);
xlabel('time');
ylabel('PDE**');

%% resample PDE
timeend = 3;
sampling = 0.001;
lengthvec = timeend/sampling+1;
for i=1:num_sim
    PDE_doub_new(:,i) = zeros(1,lengthvec);
    PDE_sing_new(:,i) = zeros(1,lengthvec);
    samp = G_R0(i) + R0(i);
    for j=1:samp
        [time_prec(:), PDE_doub_sample(:,j)] = ...
        resamplePDE(timeend, sampling, time{i,j}, PDE_doub{i,j});
        PDE_doub_new(:,i) = PDE_doub_new(:,i) + PDE_doub_sample(:,j);
        [time_prec(:), PDE_sing_sample(:,j)] = ...
            resamplePDE(timeend, sampling, time{i,j}, PDE_sing{i,j});
        PDE_sing_new(:,i) = PDE_sing_new(:,i) + PDE_sing_sample(:,j);
    end
end
time_new = linspace(0,timeend,lengthvec);

%% run backend simulation to get photocurrent (this will take a few minutes)
modelbackend = IQMmodel('Models/HSDM_backend.txtbc');
deltaJ = zeros(3001,num_sim);
time = 0:0.001:3;
for i=1:num_sim
    Etext = sprintf('interp0IQM([%s],[%s],time)',sprintf('%g,',time_new),...
    sprintf('%g,',PDE_doub_new(:,i)+0.025*PDE_sing_new(:,i)));
    % add it to the model
    indexE = variableindexIQM(modelbackend,'E');
    ms = struct(modelbackend);
    ms.variables(indexE).formula = Etext;
    modelbackend = IQMmodel(ms);
    simbackend = IQMsimulate(modelbackend,time);
    deltaJ(:,i) = simbackend.variablevalues(:,variableindexIQM(modelbackend,'deltaJ'));
end

%% plot results
figure(2);
hold on;
for i=1:num_sim
    plot(time_new, deltaJ(:,i), 'k', 'LineWidth', 1);
end
xlim([0 1]);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',15);

%% categorize events based on true number of photons
num_act = num_prec + num_nonprec;

j=1; %counter for failures
x=1; %counter for SPR
y=1; %counter for multiples

for i=1:num_sim
    if num_act(i) == 0
        truefailures(j) = i;
        j = j+1;
    elseif num_act(i) == 1
        trueSPR(x) = i;
        x = x+1;
    else
        truemult(y) = i;
        y = y+1;
    end
end

num_truefail = j-1;
num_trueSPR = x-1;
num_truemult = y-1;

%% plot true failures and their average
av_truefail = zeros(3001,1);
for i=1:num_truefail
    av_truefail = av_truefail + deltaJ(:,truefailures(i))/num_truefail;
end

figure(3);
hold on;
for i=1:num_truefail
    p1 = plot(time_new, deltaJ(:,truefailures(i)), 'k');
end
p2 = plot(time_new, av_truefail, 'r', 'LineWidth', 1.5);
ylim([-1 1]);
xlim([0 2]);
xlabel('time/s');
ylabel('\Delta J/pA');
title('True failures (66)');
legend([p1 p2], 'Stochastic', 'Average');

%% plot true SPRs and their average
av_trueSPR = zeros(3001,1);
for i=1:num_trueSPR
    av_trueSPR = av_trueSPR + deltaJ(:,trueSPR(i))/num_trueSPR;
    trueSPRs(:,i) = deltaJ(:,trueSPR(i));
end

figure(4);
hold on;
for i=1:num_trueSPR
    p1 = plot(time_new, deltaJ(:,trueSPR(i)), 'k');
end
p2 = plot(time_new, av_trueSPR, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
xlim([0 2]);
title('True SPRS (71)');
legend([p1 p2], 'Stochastic', 'Average');

max_trueSPR = max(av_trueSPR);

%% plot true mults and their average
av_truemult = zeros(3001,1);
for i=1:num_truemult
    av_truemult = av_truemult + deltaJ(:,truemult(i))/num_truemult;
end

figure(5);
hold on;
for i=1:num_truemult
    p1 = plot(time_new, deltaJ(:,truemult(i)), 'k');
end
p2 = plot(time_new, av_truemult, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
xlim([0 2]);
title('True MPRS (63)');
legend([p1 p2], 'Stochastic', 'Average');

%% calculate average
average = deltaJ(:,1)/num_sim;
for i=2:num_sim
    average = average + deltaJ(:,i)/num_sim;
end

%% coefficients of variability of ensemble
[CVarea_ALL, CVamplitude_ALL] = CV(0.001, deltaJ(1:1500,:), 0,0)

%% coefficients of variability of true SPRS
[CVarea_trueSPR, CVamplitude_trueSPR] = CV(0.001, trueSPRs(1:1500, :), 0,0)

%% find scaling factors by fitting average to each individual curve
% first compute scaling factor over rising phase
[maxR, index] = max(average); % find the peak of the response
scaling = zeros(num_sim,1);
for i=1:num_sim
    mult = average.*deltaJ(:,i);
    mult2 = average.*average;
    sum1 = sum(mult(1:index));
    sum2 = sum(mult2(1:index));
    scaling(i) = sum1/sum2;
end

%% make histogram
figure(6);
num_bins = 30;
histo = histogram(scaling, num_bins);
bins = histo.BinEdges(1:num_bins) + 0.5*histo.BinWidth;
bins = bins.';
histovalues = histo.BinCounts;
histovalues = histovalues.';
tbl = table(bins, histovalues);
xlabel('Scaling Factor');
ylabel('Counts');
set(gca,'Fontsize',20);

%% fit sum of Gaussians
% define sum of 4 Gaussians
modelfun = @(b,x) b(1)/b(2) .* exp(-((x(:,1)-b(3)).^2)./(2*b(2)^2))+...
    b(4)/b(5) .* exp(-((x(:,1)-b(6)).^2)./(2*b(5)^2))+...
    b(7)/b(8) .* exp(-((x(:,1)-b(9)).^2)./(2*b(8)^2));
% initial conditions for the parameters
beta0 = [3 0.1 0.3 ...
    1 0.2 1 ...
    0.5 0.2 2];
opts = statset('MaxIter',600);
mdl = fitnlm(tbl,modelfun,beta0, 'Options', opts)

%% plot results
figure(7);
hold on;
plot(tbl{:,1}, tbl{:,2}, 'k', 'linewidth', 1.5);
plot(tbl{:,1}, mdl.Fitted, 'r', 'linewidth', 1.5);
xlabel('Scaling factor');
ylabel('Counts');
set(gca,'Fontsize',20);
legend('Histogram', 'Sum of Gaussians Fit');

% eyeballed it: borders for SPR 0.4 to 1.4

%% filter out single photon responses
j=1; %counter for SPR
x=1; %counter for larger responses
y=1; %counter for failures
for i=1:num_sim
    if scaling(i)>0.4 
        if scaling(i)<1.4
            SPRindexes(j) = i;
            j=j+1;
        else
            larger(x) = i;
            x=x+1;
        end
    else
        smaller(y) = i;
        y = y+1;
    end
end
numSPR = j-1;
numfail = y-1;
numlarger = x-1;

%% separate responses into arrays
av_SPR = zeros(3001, 1);
for j=1:numSPR
    SPRS(:,j) = deltaJ(:,SPRindexes(j));
    av_SPR = av_SPR + deltaJ(:,SPRindexes(j))/numSPR;
end

av_fail = zeros(3001, 1);
for j=1:numfail
    failures(:,j) = deltaJ(:,smaller(j));
    av_fail = av_fail + deltaJ(:,smaller(j))/numfail;
end

av_larger = zeros(3001, 1);
for j=1:numlarger
    multiples(:,j) = deltaJ(:,larger(j));
    av_larger = av_larger + deltaJ(:,larger(j))/numlarger;
end

%% plot all SPRS
figure(8); clf;
hold on;
for i=1:numSPR
    p1 = plot(time_new, SPRS(:,i), 'k');
end
p2 = plot(time_new, av_SPR, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
xlim([0 2]);
title('Categorized SPRS (71)');
legend([p1 p2], 'Stochastic', 'Average');
            
%% plot all failures
figure(9);
hold on;
for i=1:numfail
    p1 = plot(time_new, failures(:,i), 'k');
end
p2 = plot(time_new, av_fail, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
xlim([0 2]);
title('Categorized failures (71)');
legend([p1 p2], 'Stochastic', 'Average');

%% plot multiple responses
figure(10);
hold on;
for i=1:numlarger
    p1 = plot(time_new, multiples(:,i), 'k');
end
p2 = plot(time_new, av_larger, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
xlim([0 2]);
title('Categorized MPRS (58)');
legend([p1 p2], 'Stochastic', 'Average');

%% PUBLICATION FIGURE 6
% plot all failures
figure(11); clf;
subplot(1,3,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',12);
for i=1:numfail
    p1 = plot(time_new, failures(:,i), 'k');
end
p2 = plot(time_new, av_fail, 'r', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('Photocurrent (pA)');
xlim([0 2]);
ylim([-0.5 3.5]);
title('Categorized failures');
legend([p1 p2], 'Stochastic', 'Average');

% plot all SPRS
subplot(1,3,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',12);
for i=1:numSPR
    p1 = plot(time_new, SPRS(:,i), 'k');
end
p2 = plot(time_new, av_SPR, 'r', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('Photocurrent (pA)');
xlim([0 2]);
ylim([-0.5 3.5]);
title('Categorized SPRs');

% plot multiple responses
subplot(1,3,3);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'C','Units','normalized','FontSize',12);
for i=1:numlarger
    p1 = plot(time_new, multiples(:,i), 'k');
end
p2 = plot(time_new, av_larger, 'r', 'LineWidth', 1.5);
xlabel('time (s)');
ylabel('Photocurrent (pA)');
xlim([0 2]);
ylim([-0.5 3.5]);
title('Categorized MPRs');

%% get mis-categorized responses
% look at SPRindexes and num_act
misfailindexes = [21 33 58 115 196];
misSPRindexes = [47 124 142 163 188];

misfails = zeros(3001,5);
misSPRs = zeros(3001,5);

for i=1:5
    misfails(:,i) = deltaJ(:,misfailindexes(i));
    misSPRs(:,i) = deltaJ(:,misSPRindexes(i));
end

%% PUBLICATION FIGURE S2
figure(12); clf;
subplot(1,2,1)
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',12);
for i=1:5
    plot(time_new, misfails(:,i), 'k', 'LineWidth', 1);
end
xlim([0 2]);
ylim([-0.5 3.5]);
title('Miscategorized true SPRs');
xlabel('time (s)');
ylabel('Photocurrent (pA)');
set(gca, 'FontSize', 10);

subplot(1,2,2)
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',12);
for i=1:5
    plot(time_new, misSPRs(:,i), 'k', 'LineWidth', 1);
end
xlim([0 2]);
ylim([-0.5 3.5]);
title('Miscategorized true MPRs');
xlabel('time (s)');
ylabel('Photocurrent (pA)');
set(gca, 'FontSize', 10);


%% calculate variance of failures
timestep = 0.001;
amps0 = max(failures(1:1533,:));
varamp0 = var(amps0);
area0 = sum(failures(1:1533,:))*timestep;
vararea0 = var(area0);

%% calculate CV of area and amplitude
[CVarea, CVamplitude] = CV(timestep, SPRS(1:1533,:), varamp0, vararea0)

%% calculate TTP of true SPRs and MPRs
[maxSPR, indexSPR] = max(av_trueSPR);
TTP_SPR = time_new(indexSPR)
[maxMPR, indexMPR] = max(av_truemult);
TTP_MPR = time_new(indexMPR)
