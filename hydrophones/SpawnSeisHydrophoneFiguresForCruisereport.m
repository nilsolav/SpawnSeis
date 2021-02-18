% RMS 252dB 20log10(RMS)

%% 15.02 1700
clear all
%D{1} = 'D:\DATA\S2021xxx_PH.U.SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D1_BottomMooredHydrophoneNo3\04_20210215-175947.wav';
%D{2} = 'D:\DATA\S2021xxx_PH.U.SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D2_BottomMooredHydrophoneNo4\02_20210215-180023.wav';
D{1} = 'D:\DATA\S2021xxx_PH.U.SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D1_BottomMooredHydrophoneNo3\04_20210215-182951.wav';
D{2} = 'D:\DATA\S2021xxx_PH.U.SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D2_BottomMooredHydrophoneNo4\02_20210215-182956.wav';

%% Start of pulse index
%ind = [365000 230000 724000 620000];
ind = [724000 620000];
% Plotting interval (tp) and integration interval (ti) in seconds relative 
% to the start index (ind) for the pulse
tp = [-1 2];
ti = [0 1];
cal = [.1 .1]
%% Time of start of file
for i=1:length(D)
    t0(i) = datenum(D{i}(end-18:end-4),'yyyymmdd-HHMMSS');
end

%%
close all
figure(1)
for i=1:length(D)
    % Read file
    file = D{i};
    ainf = audioinfo(file);
    dat = audioread(file);
    dt = 1/ainf.SampleRate;
    t = ((1:length(dat))-1)*dt;
    % Select pulse for plotting
    indp = (ind(i)+round(tp(1)/dt)):(ind(i)+round(tp(2)/dt));
    % Select pulse for integration
    indi = (ind(i)+round(ti(1)/dt)):(ind(i)+round(ti(2)/dt));
    % Integrate
    SEL(i) = 10*log10(trapz(t(indi),dat(indi).^2));
    
    % Plot
    subplot(1,length(D),i)
    hold on
    title(datestr(t0(i)))
    plot(t(indp)-t(ind(i)),dat(indp)*cal(i),'b')
    plot(t(indi)-t(ind(i)),dat(indi)*cal(i),'r')
    ylim([-0.01 0.0150]*cal(i))
    ylabel('Pressure (kPa)')
    xlabel('Time (s)')
end

disp(diff(SEL))

% Ca avstand mellom hydrofonar
