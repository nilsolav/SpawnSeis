%% Run SpawnSeisHydrophoneData.m to get the calibration data

%% Metata
% The data root dir depends on the computer. Place a local file called 
% 'rootdir.json' with the content {"rootdir":"d:/DATA/"} pointing to your
% local storage. Defaults to imr's file store.
clear all
if exist('rootdir.json','file')
    fid = fopen('rootdir.json','rt'); % Opening the file.
    raw = fread(fid,inf); % Reading the contents.
    fclose(fid); % Closing the file.
    str = char(raw'); % Transformation.
    par = jsondecode(str); % Using the jsondecode function to parse JSON from string.
    rootdir = par.rootdir;
else
    rootdir = '\\ces.hi.no\cruise_data';
end
% Load calibration data
load calibrationfactor.mat

%% 15.02 1700
%data\2020\S2020830_PH.U. SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D2_BottomMooredHydrophoneNo3
D{1} = fullfile(rootdir,'\2020\S2020830_PH.U. SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D1_BottomMooredHydrophoneNo2\02_20200210-163403.wav');
D{2} = fullfile(rootdir,'\2020\S2020830_PH.U. SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D2_BottomMooredHydrophoneNo3\03_20200210-163405.wav');
D{3} = fullfile(rootdir,'\2021\S2021826_PH.U.SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D14_BottomMooredHydrophone\04_20210215-182951.wav');
D{4} = fullfile(rootdir,'\2021\S2021826_PH.U.SVERDRUP II[1007]\EXPERIMENTS\HYDROPHONES\D15_BottomMooredHydrophone\02_20210215-182956.wav');
I =[1 2 14 15]; % The index in the metadata file

% Start of pulse index
ind = [585000 330000 724000 620000];
% Plotting interval (tp) and integration interval (ti) in seconds relative 
% to the start index (ind) for the pulse
tp = [-1 2];
ti = [0 1];
%cal = [.1 .1]
% Time of start of file
for i=1:length(D)
    t0(i) = datenum(D{i}(end-18:end-4),'yyyymmdd-HHMMSS');
end

%
close all
figure(1)
for i=1:length(D)
    % Read file
    file = D{i};
    ainf = audioinfo(file);
    nils = double(audioread(file,'native'));
    dat =  detrend(nils)*calibrationfactor(I(i),1);
    dt = 1/ainf.SampleRate;
    t = ((1:length(dat))-1)*dt;
    % Select pulse for plotting
    indp = (ind(i)+round(tp(1)/dt)):(ind(i)+round(tp(2)/dt));
    % Select pulse for integration
    indi = (ind(i)+round(ti(1)/dt)):(ind(i)+round(ti(2)/dt));
    % Integrate
    SEL(i) = 10*log10(trapz(t(indi),(dat(indi)*10^6).^2));% Convert Pa to uPa
    
    % Plot (should be around 100-150 Pa)
    figure(i)
    hold on
    %title(datestr(t0(i)))
    %plot(dat)
    plot(t(indp)-t(ind(i)),dat(indp),'b')
    plot(t(indi)-t(ind(i)),dat(indi),'r')
    %ylim([-50 85])
    ylabel('Pressure (Pa)')
    xlabel('Time (s)')
    F=['Fig',num2str(i),'_SEL',num2str(round(SEL(i))),'_',datestr(t0(i),30),'.png'];
    print(F,'-dpng')
end
disp(['SEL:',num2str(SEL)])
disp(diff(SEL))
%%

% Avstand frå fartøy til den første hydrofonen
R1 = 2000; % m
% Ca avstand mellom hydrofonar : 3000 m
R2 = R1 + 3000; % m

% Ca demping mellom punkta
% Es/R0  E1=Es/R1  E2=Es/R2
% E1 = Es/R1, E2=Es/R2
% E1/E2 = (Es/R1)/(Es/R2) = R2/R1
% E1/E2 = (Es/R1^2)/(Es/R2^2) = R2^2/R1^2
disp(['Actual loss between hydrophones: ',num2str(diff(SEL)),' dB re 1\muPa^2 s'])
disp(['Expected damping between hydrophones assuming spherical transmission loss: ',num2str(20*log10(R1/R2)),' dB'])
disp(['Expected damping between hydrophones assuming cylindrical transmission loss: ',num2str(10*log10(R1/R2)),' dB'])


