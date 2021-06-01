% Conversion of photovoltage into photocurrent
% Matlab funtions sin, cos work with radians by default

function [ReCurrent, ImCurrent] = conversionscript(Time, ReVoltage)

    % define parameters for conversion
    Gleak = 0.087077e-9;    % rod passive leakage conductance [Siemens]
    Cm = 13.268e-12;		% rod membrane capacitance [Farad]
    Gneg = 0.16711e-9;		% feedback conductance value at full activation [Siemens]
    tau = 0.062197;			% feedback conductance time constant [sec]
    gHat = 2.2977;          % unitless

    % make array even-sized
    arraySize = length(ReVoltage);
    if(mod(arraySize,2) ~= 0)   % uneven array size
        ReVoltage = ReVoltage(1:arraySize-1);   % cut last sample
        arraySize = arraySize-1;
    end

    % compute sample interval in frequency domain
    FSampleInterval = 1 / (Time(arraySize) - Time(1));  % caution: Matlab arrays start at 1 and end at n!

    % complex fourier transform
    VoltageFFT = fft(ReVoltage);

    % real and imaginary part
    ReVoltageFFT = real(VoltageFFT);
    ImVoltageFFT = imag(VoltageFFT); % without sign inversion

    % complex admittance
    ReAdm = zeros(arraySize, 1);
    ImAdm = zeros(arraySize, 1);
    for i=1:1:arraySize
        f = (i-1)*FSampleInterval;
        ReAdm(i) = Gleak + Gneg + Gneg * gHat / (1 + (2 * pi * f * tau)^2);
        ImAdm(i) = 2 * pi * f * (Cm - Gneg * gHat * tau / (1 + (2 * pi * f * tau)^2));
    end

    % complex product
    ReCurrentFFT = ReVoltageFFT .* ReAdm - ImVoltageFFT .* ImAdm;
    ImCurrentFFT = ReVoltageFFT .* ImAdm + ImVoltageFFT .* ReAdm;

%     load negative frequency part as expected by the complexFFT algorithm
    for i=arraySize/2 +1:1:arraySize-1
        ReCurrentFFT(i+1) = ReCurrentFFT(arraySize-i+1);
        ImCurrentFFT(i+1) = -ImCurrentFFT(arraySize-i+1);
    end

    % fourier transform back
    CurrentFFT = complex(ReCurrentFFT, ImCurrentFFT);
    Current = ifft(CurrentFFT);

    % real part: sign convention
    ReCurrent = real(Current);
    ImCurrent = imag(Current);
    ReCurrent = - ReCurrent;
end
