% Calculate coefficient of variation
% CV = sigma/mu (standard deviation/mean)
% deltaJarray must contain the separate deltaJ(t) as COLUMNS

function [CVarea, CVamplitude] = CV_new(timestep, deltaJarray, varamp0, vararea0)
   amplitudes = max(deltaJarray); % creates vector with max of each column
   meanamps = mean(amplitudes); % mean of the amplitudes
   varamps = var(amplitudes); % variance of amplitudes
   varamps = varamps - varamp0; % subtract variance of failures
   stdamps = sqrt(varamps); % standard dev
   CVamplitude = stdamps/meanamps;
   
   area = sum(deltaJarray)*timestep;
   meanarea = mean(area);
   vararea = var(area); % variance of areas
   vararea = vararea - vararea0; % subtract variance of failures
   stdarea = sqrt(vararea); % standard dev
   CVarea = stdarea/meanarea;
end