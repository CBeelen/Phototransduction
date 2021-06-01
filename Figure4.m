%% Import data
% pool data from cell 200112 and 200111
% 347 data columns (173 repeats)
% intensity 2.2 ph/mum^2
clear all; clc; close all;
fid=fopen('Exp_data/cell200112.axgt');
s1=textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'headerlines',1);
fclose(fid);
time1 = s1{1,1};
num_repeats1 = 173;

% 267 data columns (133 repeats)
% intensity 3.3 ph/mum^2
fid=fopen('Exp_data/cell200111.axgt');
s2=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
    'headerlines',1);
fclose(fid);
time2 = s2{1,1};
num_repeats2 = 133;

%% calculate mean responses
meanr1 = 1/num_repeats1*s1{1,2};
for i=3:num_repeats1+1
   meanr1 = meanr1 + 1/num_repeats1*s1{1,i}; 
end

meanr2 = 1/num_repeats2*s2{1,2};
for i=3:num_repeats2+1
   meanr2 = meanr2 + 1/num_repeats2*s2{1,i}; 
end

%% calculate scaling factors only over rising phase
% at 2660 the stimulus starts
[maxR, index] = min(meanr1); % find the peak of the response
scaling1 = zeros(num_repeats1+1,1);
for i=2:num_repeats1+1
    mult = meanr1.*s1{1,i};
    mult2 = meanr1.*meanr1;
    sum1 = sum(mult(2660:index));
    sum2 = sum(mult2(2660:index));
    scaling1(i) = sum1/sum2;
end

[maxR, index] = min(meanr2); % find the peak of the response
scaling2 = zeros(num_repeats2+1,1);
for i=2:num_repeats2+1
    mult = meanr2.*s2{1,i};
    mult2 = meanr2.*meanr2;
    sum1 = sum(mult(2660:index));
    sum2 = sum(mult2(2660:index));
    scaling2(i) = sum1/sum2;
end

%% make histogram of scaling factors for cell 200112
figure(1);
num_bins = 40;
histo = histogram(scaling1(2:num_repeats1+1), num_bins);
bins = histo.BinEdges(1:num_bins) + 0.5*histo.BinWidth;
bins = bins.';
histovalues = histo.BinCounts;
histovalues = histovalues.';
tbl = table(bins, histovalues);

%% fit sum of Gaussians for cell 200112
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

%% plot results for cell 200112
figure(2);
hold on;
plot(tbl{:,1}, tbl{:,2}, 'k', 'linewidth', 1.5);
plot(tbl{:,1}, mdl.Fitted, 'r', 'linewidth', 1.5);
xlabel('Scaling factor');
ylabel('Counts');
set(gca,'Fontsize',20);
legend('Histogram', 'Sum of Gaussians Fit');

%% find intersections for cell 200112
a1 = mdl.Coefficients{1,1};
sigma1 = mdl.Coefficients{2,1};
mu1 = mdl.Coefficients{3,1};

a2 = mdl.Coefficients{4,1};
sigma2 = mdl.Coefficients{5,1};
mu2 = mdl.Coefficients{6,1};

x = (sigma2^2*mu1 - sigma1^2*mu2 + sigma1*sigma2*sqrt((mu1-mu2)^2 ...
    + 2*(sigma2^2-sigma1^2)*log(a1*sigma2/(a2*sigma1)))) / ...
    (sigma2^2-sigma1^2)

a3 = mdl.Coefficients{7,1};
sigma3 = mdl.Coefficients{8,1};
mu3 = mdl.Coefficients{9,1};

x2 = (sigma3^2*mu2 - sigma2^2*mu3 + sigma2*sigma3*sqrt((mu2-mu3)^2 ...
    + 2*(sigma3^2-sigma2^2)*log(a2*sigma3/(a3*sigma2)))) / ...
    (sigma3^2-sigma2^2)

%% make histogram of scaling factors for cell 200111
figure(3);
num_bins = 50;
histo = histogram(scaling2(2:num_repeats2+1), num_bins);
bins = histo.BinEdges(1:num_bins) + 0.5*histo.BinWidth;
bins = bins.';
histovalues = histo.BinCounts;
histovalues = histovalues.';
tbl = table(bins, histovalues);

%% fit sum of Gaussians for cell 200111
% define sum of 4 Gaussians
modelfun = @(b,x) b(1)/b(2) .* exp(-((x(:,1)-b(3)).^2)./(2*b(2)^2))+...
    b(4)/b(5) .* exp(-((x(:,1)-b(6)).^2)./(2*b(5)^2))+...
    b(7)/b(8) .* exp(-((x(:,1)-b(9)).^2)./(2*b(8)^2))+...
    b(10)/b(11) .* exp(-((x(:,1)-b(12)).^2)./(2*b(11)^2));
% initial conditions for the parameters
beta0 = [3.3 0.4 0 ...
    0.7 0.1 1 ...
    0.7 0.1 1.5 ...
    1 0.1 2];
mdl = fitnlm(tbl,modelfun,beta0)

%% plot results for cell 200111
figure(4);
hold on;
plot(tbl{:,1}, tbl{:,2}, 'k', 'linewidth', 1.5);
plot(tbl{:,1}, mdl.Fitted, 'r', 'linewidth', 1.5);
xlabel('Scaling factor');
ylabel('Counts');
set(gca,'Fontsize',20);
legend('Histogram', 'Sum of Gaussians Fit');

%% intersections for 200111: mixture of intersections and "manual" estimate

%% filter out single photon responses
% for cell 200112: SPRs between 0.12 and 0.93
% for cell 200111: SPRs between 0.84 and 1.38
j=1; %counter for SPR
x=1; %counter for larger responses
y=1; %counter for failures
for i=2:num_repeats1+1
    if scaling1(i)>0.12 
        if scaling1(i)<0.93
            SPRindexes1(j) = i;
            j=j+1;
        else
            larger1(x) = i;
            x=x+1;
        end
    else
        smaller1(y) = i;
        y = y+1;
    end
end
numSPR1 = j-1;
numfail1 = y-1;
numlarger1 = x-1;

% reset counters
j=1; %counter for SPR
x=1; %counter for larger responses
y=1; %counter for failures
for i=2:num_repeats2+1
    if scaling2(i)>0.84 
        if scaling2(i)<1.38
            SPRindexes2(j) = i;
            j=j+1;
        else
            larger2(x) = i;
            x=x+1;
        end
    else
        smaller2(y) = i;
        y = y+1;
    end
end
numSPR2 = j-1;
numfail2 = y-1;
numlarger2 = x-1;

%% separate responses into arrays
for j=1:numSPR1
    SPRS1(:,j) = s1{1,SPRindexes1(j)};
end

numneg1 = 0;
for j=1:numfail1
    failures1(:,j) = s1{1,smaller1(j)};
    if scaling1(smaller1(j)) < 0
        negatives1(:,j) = s1{1,smaller1(j)};
        numneg1 = numneg1 +1;
    end
end

for j=1:numlarger1
    multiples1(:,j) = s1{1,larger1(j)};
end

for j=1:numSPR2
    SPRS2(:,j) = s2{1,SPRindexes2(j)};
end

numneg2 = 0;
for j=1:numfail2
    failures2(:,j) = s2{1,smaller2(j)};
    if scaling2(smaller2(j)) < 0
        negatives2(:,j) = s2{1,smaller2(j)};
        numneg2 = numneg2 +1;
    end
end

for j=1:numlarger2
    multiples2(:,j) = s2{1,larger2(j)};
end

%% convert to photocurrent
for i=1:numSPR1
    [ReCurrentSPR1(:,i), ImCurrentSPR1(:,i)] = conversionscript(s1{1,1}, SPRS1(:,i));
end
for i=1:numfail1
    [ReCurrentfail1(:,i), ImCurrentfail1(:,i)] = conversionscript(s1{1,1}, failures1(:,i));
end
for i=1:numlarger1
    [ReCurrentmult1(:,i), ImCurrentmult1(:,i)] = conversionscript(s1{1,1}, multiples1(:,i));
end

for i=1:numSPR2
    [ReCurrentSPR2(:,i), ImCurrentSPR2(:,i)] = conversionscript(s2{1,1}, SPRS2(:,i));
end
for i=1:numfail2
    [ReCurrentfail2(:,i), ImCurrentfail2(:,i)] = conversionscript(s2{1,1}, failures2(:,i));
end
for i=1:numlarger2
    [ReCurrentmult2(:,i), ImCurrentmult2(:,i)] = conversionscript(s2{1,1}, multiples2(:,i));
end

%% filter

for i=1:numSPR1
    ReCurrentSPR1(:,i) = lowpass(ReCurrentSPR1(:,i),40,5000);
end
for i=1:numfail1
    ReCurrentfail1(:,i) = lowpass(ReCurrentfail1(:,i),40,5000);
end
for i=1:numlarger1
    ReCurrentmult1(:,i) = lowpass(ReCurrentmult1(:,i),40,5000);
end

for i=1:numSPR2
    ReCurrentSPR2(:,i) = lowpass(ReCurrentSPR2(:,i),40,5000);
end
for i=1:numfail2
    ReCurrentfail2(:,i) = lowpass(ReCurrentfail2(:,i),40,5000);
end
for i=1:numlarger2
    ReCurrentmult2(:,i) = lowpass(ReCurrentmult2(:,i),40,5000);
end

%% average single photon responses etc
avSPR1 = zeros(length(s1{1,1}),1);
for i=1:numSPR1
    avSPR1 = avSPR1 + ReCurrentSPR1(:,i)./numSPR1;
end

avfail1 = zeros(length(s1{1,1}),1);
for i=1:numfail1
    avfail1 = avfail1 + ReCurrentfail1(:,i)./numfail1;
end

avmult1 = zeros(length(s1{1,1}),1);
for i=1:numlarger1
    if i ~= 15      % i= 15 is an outlier
        avmult1 = avmult1 + ReCurrentmult1(:,i)./(numlarger1-1);
    end
end

avSPR2 = zeros(length(s2{1,1}),1);
for i=1:numSPR2
    avSPR2 = avSPR2 + ReCurrentSPR2(:,i)./numSPR2;
end

avfail2 = zeros(length(s2{1,1}),1);
for i=1:numfail2
    avfail2 = avfail2 + ReCurrentfail2(:,i)./numfail2;
end

avmult2 = zeros(length(s2{1,1}),1);
for i=1:numlarger2
    avmult2 = avmult2 + ReCurrentmult2(:,i)./numlarger2;
end

%% find SPR amps
SPRamp1 = max(avSPR1(2660:5000));
SPRamp2 = max(avSPR2(2660:5000));

%% plot average SPRs
figure(5);
hold on;
plot(s1{1,1}, avSPR1/SPRamp1, 'k');
plot(s2{1,1}, avSPR2/SPRamp2, 'r');
xlabel('time/s');
ylabel('scaled \Delta J');
set(gca,'Fontsize',20);
xlim([0.3,1.5]);
title('Average SPR');

%% plot average failures
figure(6);
hold on;
plot(s1{1,1}, avfail1/SPRamp1, 'k');
plot(s2{1,1}, avfail2/SPRamp2, 'r');
xlabel('time/s');
ylabel('scaled \Delta J');
set(gca,'Fontsize',20);
xlim([0.3,1.5]);
title('Average failures');

%% plot all SPRs
figure(7);
hold on;
for i=1:numSPR1
    plot(s1{1,1}, ReCurrentSPR1(:,i)/SPRamp1, 'k');
end
for i=1:numSPR2
    plot(s2{1,1}, ReCurrentSPR2(:,i)/SPRamp2, 'r');
end
xlabel('time/s');
ylabel('scaled \Delta J');
set(gca,'Fontsize',20);
xlim([0.3,1.5]);
title('All SPRs');

%% paste together
ReCurrentSPR = [ReCurrentSPR1/SPRamp1, ReCurrentSPR2/SPRamp2];
ReCurrentfail = [ReCurrentfail1/SPRamp1, ReCurrentfail2/SPRamp2];
ReCurrentmult = [ReCurrentmult1/SPRamp1, ReCurrentmult2/SPRamp2];
numSPR = numSPR1 + numSPR2;
numfail = numfail1 + numfail2;
numlarger = numlarger1+numlarger2;

%% new averages
avSPR = zeros(length(s1{1,1}),1);
for i=1:numSPR
    avSPR = avSPR + ReCurrentSPR(:,i)./numSPR;
end

avfail = zeros(length(s1{1,1}),1);
for i=1:numfail
    avfail = avfail + ReCurrentfail(:,i)./numfail;
end

avmult = zeros(length(s1{1,1}),1);
for i=1:numlarger
    if i ~= 15      % i= 15 is an outlier
        avmult = avmult + ReCurrentmult(:,i)./(numlarger-1);
    end
end

%% calculate CV
timestep = s1{1,1}(2)-s1{1,1}(1);
% calculate variance of failures
amps0 = min(ReCurrentfail(2663:10000,:));
varamp0 = var(amps0);
area0 = sum(ReCurrentfail(2663:10000,:))*timestep;
vararea0 = var(area0);
[CVarea, CVamplitude] = CV(timestep, ReCurrentSPR(2660:10000,:), varamp0, vararea0)

%% TTP
[maxR, index] = max(avSPR(1000:5000));
time_LED = 0.5326;
TTP = s1{1,1}(index+1000)-time_LED


%% plot categorized responses
figure(8);clf;
subplot(1,3,1);
hold on;
pbaspect([1 1 1]);
xlim([-0.2,1.3]);
ylim([-2.5 3]);
for i=1:numfail
    p1 = plot(s1{1,1}-0.53, ReCurrentfail(:,i), 'k');
end
p2 = plot(s1{1,1}-0.53, avfail, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('Scaled \Delta J');
%set(gca,'Fontsize',20);
legend([p1 p2], {'All traces', 'Average'});
title('Failures');

subplot(1,3,2);
hold on;
pbaspect([1 1 1]);
xlim([-0.2,1.3]);
ylim([-1.5 3]);
for i=1:numSPR
    p1 = plot(s1{1,1}-0.53, ReCurrentSPR(:,i), 'k');
end
p2 = plot(s1{1,1}-0.53, avSPR, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('Scaled \Delta J');
%set(gca,'Fontsize',20);
%legend([p1 p2], {'SPR traces', 'average SPR trace'});
title('SPRs');

subplot(1,3,3);
hold on;
pbaspect([1 1 1]);
xlim([-0.2,1.3]);
ylim([-1.5 6.5]);
for i=1:numlarger
    if i ~= 15      % i= 15 is an outlier
        p1 = plot(s1{1,1}-0.53, ReCurrentmult(:,i), 'k');
    end
end
p2 = plot(s1{1,1}-0.53, avmult, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('Scaled \Delta J');
title('MPRs');
%set(gca,'Fontsize',20);
%legend([p1 p2], {'MPR traces', 'average MPR trace'});

%% Figure 4 in the publication
figure(9);clf;
subplot(1,3,1);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'A','Units','normalized','FontSize',15);
xlim([-0.2,1.3]);
ylim([-1.5 6.5]);
for i=1:numfail
    p1 = plot(s1{1,1}-0.5326, ReCurrentfail(:,i), 'k');
end
p2 = plot(s1{1,1}-0.5326, avfail, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('Scaled {\Delta}J');
set(gca,'Fontsize',15);
legend([p1 p2], {'All traces', 'Average'});
title('Failures');

subplot(1,3,2);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'B','Units','normalized','FontSize',15);
xlim([-0.2,1.3]);
ylim([-1.5 6.5]);
for i=1:numSPR
    p1 = plot(s1{1,1}-0.5326, ReCurrentSPR(:,i), 'k');
end
p2 = plot(s1{1,1}-0.5326, avSPR, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('Scaled {\Delta}J');
set(gca,'Fontsize',15);
%legend([p1 p2], {'SPR traces', 'average SPR trace'});
title('SPRs');

subplot(1,3,3);
hold on;
pbaspect([1 1 1]);
text(0.025,0.95,'C','Units','normalized','FontSize',15);
xlim([-0.2,1.3]);
ylim([-1.5 6.5]);
for i=1:numlarger
    if i ~= 15      % i= 15 is an outlier
        p1 = plot(s1{1,1}-0.5326, ReCurrentmult(:,i), 'k');
    end
end
p2 = plot(s1{1,1}-0.5326, avmult, 'r', 'LineWidth', 1.5);
xlabel('time/s');
ylabel('Scaled {\Delta}J');
title('MPRs');
set(gca,'Fontsize',15);
%legend([p1 p2], {'MPR traces', 'average MPR trace'});

%% save average SPR trace
save('Exp_data/avSPR_exp.mat', 'time1', 'avSPR');