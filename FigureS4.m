%% Import data
% 347 data columns (173 repeats)
% intensity 2.2 ph/mum^2
clear all; clc; close all;
fid=fopen('Exp_data/cell200112.axgt');
s=textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'headerlines',1);
fclose(fid);
time = s{1,1};
num_repeats = 173;

%% calculate mean response
meanr = 1/num_repeats*s{1,2};
for i=3:num_repeats+1
   meanr = meanr + 1/num_repeats*s{1,i}; 
end

%% calculate scaling factor only over rising phase
% at 2660 the stimulus starts
[maxR, index] = min(meanr); % find the peak of the response
scaling3 = zeros(num_repeats+1,1);
for i=2:num_repeats+1
    mult = meanr.*s{1,i};
    mult2 = meanr.*meanr;
    sum1 = sum(mult(2660:index));
    sum2 = sum(mult2(2660:index));
    scaling3(i) = sum1/sum2;
end

%% make histogram
figure(1);clf;
num_bins = 40;
histo = histogram(scaling3(2:num_repeats+1), num_bins);
bins = histo.BinEdges(1:num_bins) + 0.5*histo.BinWidth;
bins = bins.';
histovalues = histo.BinCounts;
histovalues = histovalues.';
tbl = table(bins, histovalues);
xlabel('Scaling Factor');
ylabel('Counts');
set(gca,'Fontsize',15);
xlim([-1 4]);
ylim([0 15]);

%% fit sum of Gaussians
% define sum of 4 Gaussians
modelfun = @(b,x) b(1)/b(2) .* exp(-((x(:,1)-b(3)).^2)./(2*b(2)^2))+...
    b(4)/b(5) .* exp(-((x(:,1)-b(6)).^2)./(2*b(5)^2))+...
    b(7)/b(8) .* exp(-((x(:,1)-b(9)).^2)./(2*b(8)^2))+...
    b(10)/b(11) .* exp(-((x(:,1)-b(12)).^2)./(2*b(11)^2))+...
    b(13)/b(14) .* exp(-((x(:,1)-b(15)).^2)./(2*b(14)^2));
% initial conditions for the parameters
beta0 = [0.5 0.1 0 ...
    3 0.2 0.5 ...
    0.5 0.1 1 ...
    1 0.1 1.5 ...
    0.5 0.1 2];
opts = statset('MaxIter',600);
mdl = fitnlm(tbl,modelfun,beta0, 'Options', opts)

%%
time_fit = -1:0.01:4;
exp_fit = mdl.Coefficients{1,1}/mdl.Coefficients{2,1}.* exp(-((time_fit-mdl.Coefficients{3,1}).^2)./(2*mdl.Coefficients{2,1}^2))+...
    mdl.Coefficients{4,1}/mdl.Coefficients{5,1} .* exp(-((time_fit-mdl.Coefficients{6,1}).^2)./(2*mdl.Coefficients{5,1}^2))+...
    mdl.Coefficients{7,1}/mdl.Coefficients{8,1} .* exp(-((time_fit-mdl.Coefficients{9,1}).^2)./(2*mdl.Coefficients{8,1}^2))+...
    mdl.Coefficients{10,1}/mdl.Coefficients{11,1} .* exp(-((time_fit-mdl.Coefficients{12,1}).^2)./(2*mdl.Coefficients{11,1}^2))+...
    mdl.Coefficients{13,1}/mdl.Coefficients{14,1} .* exp(-((time_fit-mdl.Coefficients{15,1}).^2)./(2*mdl.Coefficients{14,1}^2));


%% plot all in one for publication
figure(2); clf;
subplot(2,3,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',12);
plot(time-0.53, -s{2}/1e-3, 'k');
plot(time-0.53, -s{3}/1e-3, 'k');
plot(time-0.53, -s{4}/1e-3, 'k');
plot(time-0.53, -s{5}/1e-3, 'k');
plot(time-0.53, -s{6}/1e-3, 'k');
plot(time-0.53, -s{7}/1e-3, 'k');
plot(time-0.53, -s{8}/1e-3, 'k');
plot(time-0.53, -s{9}/1e-3, 'k');
plot(time-0.53, -s{10}/1e-3, 'k');
plot(time-0.53, -s{11}/1e-3, 'k');
plot(time-0.53, -s{12}/1e-3, 'k');
plot(time-0.53, -s{13}/1e-3, 'k');
plot(time-0.53, -s{14}/1e-3, 'k');
plot(time-0.53, -s{15}/1e-3, 'k');
plot(time-0.53, -s{16}/1e-3, 'k');
plot(time-0.53, -s{17}/1e-3, 'k');
plot(time-0.53, -s{18}/1e-3, 'k');
p1 = plot(time-0.53, -s{19}/1e-3, 'k');
p2 = plot(time-0.53, -meanr/1e-3, 'r', 'LineWidth', 1.2);
xlabel('time (s)');
ylabel('Photovoltage (mV)');
xlim([-0.2 1.5]);
ylim([-1.2 4.5]);
set(gca,'Fontsize',10);
legend([p1 p2], 'Recordings', 'Average');
title('Compute average response');

subplot(2,3,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',12);
i=4; % 10, 14 and 4 look nice
xlim([-0.1,1.5]);
ylim([-1,4.5]);
plot(time-0.53, -s{1,i}/1e-3, 'k', 'lineWidth', 1.5);
plot(time-0.53, -scaling3(i)*meanr/1e-3, 'r', 'lineWidth', 1.5);
xlabel('time (s)');
ylabel('Photovoltage (mV)');
set(gca,'Fontsize',10);
legend('Recording', 'Fitted Mean');
title('Scale mean...')

subplot(2,3,3);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'C','Units','normalized','FontSize',12);
i=14; % 10, 14 and 4 look nice
xlim([-0.1,1.5]);
ylim([-1,4.5]);
plot(time-0.53, -s{1,i}/1e-3, 'k', 'lineWidth', 1.5);
plot(time-0.53, -scaling3(i)*meanr/1e-3, 'r', 'lineWidth', 1.5);
xlabel('time (s)');
ylabel('Photovoltage (mV)');
set(gca,'Fontsize',10);
legend('Recording', 'Fitted Mean');
title('...to each...')

subplot(2,3,4);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'D','Units','normalized','FontSize',12);
i=10; % 10, 14 and 4 look nice
xlim([-0.1,1.5]);
ylim([-1,4.5]);
plot(time-0.53, -s{1,i}/1e-3, 'k', 'lineWidth', 1.5);
plot(time-0.53, -scaling3(i)*meanr/1e-3, 'r', 'lineWidth', 1.5);
xlabel('time (s)');
ylabel('Photovoltage (mV)');
set(gca,'Fontsize',10);
legend('Recording', 'Fitted Mean');
title('...recording');

subplot(2,3,5);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'E','Units','normalized','FontSize',12);
num_bins = 40;
histo = histogram(scaling3(2:num_repeats+1), num_bins);
bins = histo.BinEdges(1:num_bins) + 0.5*histo.BinWidth;
bins = bins.';
histovalues = histo.BinCounts;
histovalues = histovalues.';
tbl = table(bins, histovalues);
xlabel('Scaling Factor');
ylabel('Counts');
set(gca,'Fontsize',10);
xlim([-1 4]);
ylim([0 15]);
title('Histogram of scaling factors');

subplot(2,3,6)
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'F','Units','normalized','FontSize',12);
text(0.02,0.6,'Failures','Units','normalized','FontSize',10);
text(0.24,0.97,'SPRs','Units','normalized','FontSize',10);
text(0.6,0.7,'MPRs','Units','normalized','FontSize',10);
plot(tbl{:,1}, tbl{:,2}, 'k', 'linewidth', 1.5);
plot(time_fit, exp_fit, 'r', 'linewidth', 1.5);
xline(0.12, '--', 'LineWidth', 1.5)
xline(0.93, '--', 'LineWidth', 1.5)
xlabel('Scaling factor');
ylabel('Counts');
set(gca,'Fontsize',10);
legend('Histogram', 'Fit', 'Intersections');
xlim([-1 4]);
ylim([0 15]);
title('Fit Sum of Gaussians to Histogram');

%% filter out single photon responses
j=1; %counter for SPR
x=1; %counter for larger responses
y=1; %counter for failures
for i=2:num_repeats+1
    if scaling3(i)>0.12 
        if scaling3(i)<0.93
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
for j=1:numSPR
    SPRS(:,j) = s{1,SPRindexes(j)};
end

numneg = 0;
for j=1:numfail
    failures(:,j) = s{1,smaller(j)};
    if scaling3(smaller(j)) < 0
        negatives(:,j) = s{1,smaller(j)};
        numneg = numneg +1;
    end
end

for j=1:numlarger
    multiples(:,j) = s{1,larger(j)};
end
%% plot all SPRS
figure(9);
hold on;
for i=1:numSPR
plot(s{1,1}, SPRS(:,i));
end
            
%% plot all failures
figure(10);
hold on;
for i=1:numfail
plot(s{1,1}, failures(:,i));
end

%% plot multiple responses
figure(11);
hold on;
for i=1:numlarger
plot(s{1,1}, multiples(:,i));
end

%% plot negative scaling factor
figure(12);
hold on;
for i=1:numneg
plot(s{1,1}, negatives(:,i));
end

%% plot all in one
figure(13);
hold on;
for i=1:numSPR
plot(s{1,1}, SPRS(:,i), 'r');
end
hold on;
for i=1:numfail
plot(s{1,1}, failures(:,i), 'k');
end
for i=1:numlarger
plot(s{1,1}, multiples(:,i), 'g');
end

%% define range for CV: stimulus start to end of recording
timestep = s{1,1}(2)-s{1,1}(1);
%% calculate variance of failures
amps0 = min(failures(2660:end,:));
varamp0 = var(amps0);
area0 = sum(failures(2660:end,:))*timestep;
vararea0 = var(area0);
%% calculate CV of area and amplitude
[CVarea, CVamplitude] = CV_new(timestep, SPRS(2660:end,:), varamp0, vararea0)

%% nice plot of SPRS etc
figure(14);
hold on;
for i=1:numSPR
plot(s{1,1}, SPRS(:,i), 'r');
end
hold on;
for i=1:numfail
plot(s{1,1}, failures(:,i), 'k');
end
for i=1:numlarger
plot(s{1,1}, multiples(:,i), 'g');
end
xlabel('time/s');
ylabel('Photovoltage');
set(gca,'Fontsize',20);
xlim([0.5,1.5]);

%% average single photon response
avSPRV = zeros(length(s{1,1}),1);
for i=1:numSPR
    avSPRV = avSPRV + SPRS(:,i)./numSPR;
end

figure(15);
hold on;
plot(s{1,1}, avSPRV);
xlabel('time/s');
ylabel('Photovoltage');
set(gca,'Fontsize',20);
xlim([0.5 1.75]);

%% convert to photocurrent
for i=1:numSPR
    [ReCurrentSPR(:,i), ImCurrentSPR(:,i)] = conversionscript(s{1,1}, SPRS(:,i));
end
for i=1:numfail
    [ReCurrentfail(:,i), ImCurrentfail(:,i)] = conversionscript(s{1,1}, failures(:,i));
end
for i=1:numlarger
    [ReCurrentmult(:,i), ImCurrentmult(:,i)] = conversionscript(s{1,1}, multiples(:,i));
end

% %% plot again
% figure(16);
% hold on;
% for i=1:numSPR
% plot(s{1,1}, ReCurrentSPR(:,i)*1e12, 'k');
% end
% hold on;
% for i=1:numfail
% plot(s{1,1}, ReCurrentfail(:,i)*1e12, 'k');
% end
% for i=1:numlarger
% plot(s{1,1}, ReCurrentmult(:,i)*1e12, 'k');
% end
% xlabel('time/s');
% ylabel('\Delta J/pA');
% set(gca,'Fontsize',20);
% xlim([0.5,1.5]);
% ylim([-1 2]);

%% average single photon response
avSPR = zeros(length(s{1,1}),1);
for i=1:numSPR
    avSPR = avSPR + ReCurrentSPR(:,i)./numSPR;
end

avfail = zeros(length(s{1,1}),1);
for i=1:numfail
    avfail = avfail + ReCurrentfail(:,i)./numfail;
end

avmult = zeros(length(s{1,1}),1);
for i=1:numlarger
    if i ~= 15      % i= 15 is an outlier
        avmult = avmult + ReCurrentmult(:,i)./(numlarger-1);
    end
end


figure(17);
hold on;
plot(s{1,1}, avSPR);
xlabel('time/s');
ylabel('\Delta J/A');
set(gca,'Fontsize',20);
xlim([0.5,1.5]);
ylim([-0.1e-12,1e-12]);
title('Average SPR');

%% compare to simulation
modeloverall = IQMmodel('../../Beelen_new/models/Beelen_new_SPR.txtbc');
time = 0:0.001:1.5;
simoverall = IQMPsimulate(modeloverall,time);
deltaJoverall = simoverall.variablevalues(:,variableindexIQM(modeloverall,'deltaJ'));

% %% plot comparison
% figure(18);
% hold on;
% plot(s{1,1}, avSPR/max(avSPR(2660:8000)), 'r', 'LineWidth', 1.5);
% plot(time+0.53, deltaJoverall/max(deltaJoverall), 'k', 'LineWidth', 1.5);
% legend('Measured', 'Simulation');
% xlim([0.53 1.75]);
% ylim([-0.1 1.1]);
% xlabel('time/s','Fontsize', 15);
% ylabel('Photocurrent (normalized)','Fontsize', 15);
% set(gca,'Fontsize',20);


% %% plot comparison
% figure(19);
% hold on;
% plot(s{1,1}, avSPR, 'r');
% plot(time+0.53, deltaJoverall*1e-12, 'k');
% legend('Photocurrent', 'Simulation');
% xlim([0.5 1.75]);
% ylim([-0.1e-12 1.1e-12]);
% xlabel('time/s','Fontsize', 15);
% ylabel('Photocurrent/A','Fontsize', 15);
% set(gca,'Fontsize',10);

% %% comparison 2
% figure(20);
% hold on;
% plot(s{1,1}, avSPR/max(avSPR(2660:8000)), 'r');
% plot(time+0.53, deltaJoverall/max(deltaJoverall), 'k');
% plot(s{1,1}, avSPRV/min(avSPRV(2660:8000)), 'g');
% legend('Photocurrent', 'Simulation', 'Photovoltage');
% xlim([0.5 1.75]);
% ylim([-0.1 1.1]);
% xlabel('time/s','Fontsize', 15);
% ylabel('Photocurrent (normalized)','Fontsize', 15);
% set(gca,'Fontsize',10);


%% define range for CV: stimulus start to end of recording
timestep = s{1,1}(2)-s{1,1}(1);
%% calculate variance of failures
amps0 = min(ReCurrentfail(2660:5000,:));
varamp0 = var(amps0);
area0 = sum(ReCurrentfail(2660:5000,:))*timestep;
vararea0 = var(area0);
%% calculate CV of area and amplitude
[CVarea, CVamplitude] = CV_new(timestep, ReCurrentSPR(2660:5000,:), varamp0, vararea0)

%%
 amplitudes = max(ReCurrentSPR(2660:5000,:)); % creates vector with max of each column
   meanamps = mean(amplitudes); % mean of the amplitudes
   varamps = var(amplitudes); % variance of amplitudes
   varamps2 = varamps - varamp0; % subtract variance of failures
   stdamps = sqrt(varamps2); % standard dev
   CVamplitude = stdamps/meanamps;
   
   area = sum(ReCurrentSPR(2660:5000,:))*timestep;
   meanarea = mean(area);
   vararea = var(area); % variance of areas
   vararea2 = vararea - vararea0; % subtract variance of failures
   stdarea = sqrt(vararea2); % standard dev
   CVarea = stdarea/meanarea;
% %% plot all SPRS
% figure(21);
% hold on;
% for i=1:numSPR
% plot(s{1,1}, ReCurrentSPR(:,i));
% end
% xlim([0.5,1.5]);
% ylim([-1e-12,3e-12]);


%%
% figure(21);
% hold on;
% for i=175:340
%     plot(s{1,1}, s{1,i});
% end
% 
% for i=2:174
%     plot(s{1,1}, -4*s{1,i});
% end
% 
% plot(s{1,1}, avSPR*0.5e11);
% 
% %%
% figure(22);
% hold on;
% plot(s{1,1}, s{1,175});
% for i=1:numSPR
%     plot(s{1,1}, ReCurrentSPR(:,i)*1e10, 'k');
% end
% plot(s{1,1}, avSPR*1e10, 'r');
% plot(time+0.53, deltaJoverall*1e-2, 'g');

%% plot for publication Fig2
figure(25);
hold on;
for i=1:numSPR
    plot(s{1,1}, ReCurrentSPR(:,i)*1e12, 'r');
end
for i=1:numfail
    plot(s{1,1}, ReCurrentfail(:,i)*1e12, 'k');
end
for i=1:numlarger
    plot(s{1,1}, ReCurrentmult(:,i)*1e12, 'g');
end
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
xlim([0.53,1.5]);
ylim([-1 3]);
hold off;

% create smaller axes in top right, and plot on it
axes('Position',[0.6 0.6 .3 .3]);
box on;
hold on;
xlim([0.2,1.8]);
ylim([-2 7]);
for i=1:numSPR
    plot(s{1,1}, ReCurrentSPR(:,i)*1e12, 'k');
end
for i=1:numfail
    plot(s{1,1}, ReCurrentfail(:,i)*1e12, 'k');
end
for i=1:numlarger
    plot(s{1,1}, ReCurrentmult(:,i)*1e12, 'k');
end
hold off;

%% new plots for publication
figure(26);
hold on;
xlim([-0.2,1.3]);
ylim([-0.75 1.75]);
for i=1:numSPR
    p1 = plot(s{1,1}-0.53, ReCurrentSPR(:,i)*1e12, 'k');
end
p2 = plot(s{1,1}-0.53, avSPR*1e12, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
legend([p1 p2], {'SPR traces', 'average SPR trace'});


%%
figure(27);
hold on;
xlim([-0.2,1.3]);
ylim([-1 1.5]);
for i=1:numfail
    p1 = plot(s{1,1}-0.53, ReCurrentfail(:,i)*1e12, 'k');
end
p2 = plot(s{1,1}-0.53, avfail*1e12, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
legend([p1 p2], {'Failure traces', 'average failure trace'});

%%
figure(28);
hold on;
xlim([-0.2,1.3]);
ylim([-1 3.5]);
for i=1:numlarger
    if i ~= 15      % i= 15 is an outlier
        p1 = plot(s{1,1}-0.53, ReCurrentmult(:,i)*1e12, 'k');
    end
end
p2 = plot(s{1,1}-0.53, avmult*1e12, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
legend([p1 p2], {'MPR traces', 'average MPR trace'});

%% low pass filter
ReCurrentSPR_filt = zeros(length(s{1,1}), numSPR);
for i=1:numSPR
    ReCurrentSPR_filt(:,i) = lowpass(ReCurrentSPR(:,i),40,5000);
end
ReCurrentfail_filt = zeros(length(s{1,1}), numfail);
for i=1:numfail
    ReCurrentfail_filt(:,i) = lowpass(ReCurrentfail(:,i),40,5000);
end
ReCurrentmult_filt = zeros(length(s{1,1}), numlarger-1);
for i=1:numlarger
    ReCurrentmult_filt(:,i) = lowpass(ReCurrentmult(:,i),40,5000);
end

%% new averages
avSPR_filt = zeros(length(s{1,1}),1);
for i=1:numSPR
    avSPR_filt = avSPR_filt + ReCurrentSPR_filt(:,i)./numSPR;
end
avfail_filt = zeros(length(s{1,1}),1);
for i=1:numfail
    avfail_filt = avfail_filt + ReCurrentfail_filt(:,i)./numfail;
end
avmult_filt = zeros(length(s{1,1}),1);
for i=1:numlarger
    if i ~= 15  % i=15 is an outlier
        avmult_filt = avmult_filt + ReCurrentmult_filt(:,i)./(numlarger-1);
    end
end

%% calculate variance of failures
amps0 = min(ReCurrentfail_filt(2660:5000,:));
varamp0 = var(amps0);
area0 = sum(ReCurrentfail_filt(2660:5000,:))*timestep;
vararea0 = var(area0);
%% calculate CV of area and amplitude
[CVarea, CVamplitude] = CV_new(timestep, ReCurrentSPR_filt(2660:5000,:), varamp0, vararea0)

%% filtered
figure(30);
hold on;
xlim([-0.2,1.3]);
ylim([-0.75 1.75]);
for i=1:numSPR
    p1 = plot(s{1,1}-0.53, ReCurrentSPR_filt(:,i)*1e12, 'k');
end
p2 = plot(s{1,1}-0.53, avSPR_filt*1e12, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
legend([p1 p2], {'SPR traces', 'average SPR trace'});

%%
figure(31);
hold on;
xlim([-0.2,1.3]);
ylim([-1 1.5]);
for i=1:numfail
    p1 = plot(s{1,1}-0.53, ReCurrentfail_filt(:,i)*1e12, 'k');
end
p2 = plot(s{1,1}-0.53, avfail_filt*1e12, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
legend([p1 p2], {'Failure traces', 'average failure trace'});

%%
figure(32);
hold on;
xlim([-0.2,1.3]);
ylim([-1 3.5]);
for i=1:numlarger
    if i ~= 15      % i= 15 is an outlier
        p1 = plot(s{1,1}-0.53, ReCurrentmult_filt(:,i)*1e12, 'k');
    end
end
p2 = plot(s{1,1}-0.53, avmult_filt*1e12, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('\Delta J/pA');
set(gca,'Fontsize',20);
legend([p1 p2], {'MPR traces', 'average MPR trace'});

