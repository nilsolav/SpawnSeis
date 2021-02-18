% Hydrophone measurements for the SpawnSeis exposure project
% Requires the m_map package

%% Metata
% The data root dir depends on the computer. Place a local file called 
% 'rootdir.json' with the content {"rootdir":"d:/DATA/"} pointing to your
% local storage. Defaults to imr's file store.

if exist('rootdir.json','file')
    fid = fopen('rootdir.json','rt'); % Opening the file.
    raw = fread(fid,inf); % Reading the contents.
    fclose(fid); % Closing the file.
    str = char(raw'); % Transformation.
    par = jsondecode(str); % Using the jsondecode function to parse JSON from string.
    rootdir = par.rootdir;
else
    rootdir = '\\ces.hi.no\';
end

% Metadata for each hydrophone 
[~,~,Hmeta_raw] = xlsread('SpawnSeisHydrophoneMetaData.csv');
Hmeta=cell2struct(Hmeta_raw(2:end,[1:7]),Hmeta_raw(1,[1:7]),2);

% Metadata for each hydrophone deployment
[~,~,Dmeta_raw] = xlsread('SpawnSeisHydrophoneDeploymentMetaData.csv');
Dmeta=cell2struct(Dmeta_raw(2:end,:),Dmeta_raw(1,:),2);

% Get metadata for the treatments
[~,~,Tmeta_raw] = xlsread('treatments.csv');
Tmeta=cell2struct(Tmeta_raw(2:end,:),Tmeta_raw(1,:),2);

%% Get calibration data per deployment
for i=[1 2 5 6]% 6 7 8] %1:length(Dmeta)
   knownPa = 10^((Dmeta(i).CalibrationLevel/20)-6); % Pa reference pressure rms-value in calibrator with naxys coupler (Se doc in H2... folder)
   % Beginning calibration
   calfile = Dmeta(i).CalibrationFileBeginning; %nok med ei fil
   dat_temp =audioread(calfile,'native'); %les inn rå-versjonen, utan native blir normalisert til 1.
   dat=detrend(double(dat_temp));
%   figure(i)
%   plot(dat)
%   hold on;
%   title(Dmeta(i).Folder, 'Interpreter', 'none')

   
%    if strcmp(Dmeta(i).CalibrationStopIndBeginning,'end')
%        dat = dat_temp(Dmeta(i).CalibrationStartIndBeginning:end);
%    else
%        dat = dat_temp(Dmeta(i).CalibrationStartIndBeginning:Dmeta(i).CalibrationStopIndBeginning);       
%    end
   
   dat_inf = audioinfo(calfile);
   caldata(i,1) = knownPa/rms(dat);
   disp(['This value shoud be 26.6. Check: ',num2str(rms(dat)*caldata(i,1))])
   Fs=dat_inf.SampleRate;

   
   rmsValueCal=rms(dat);
   calibrationfactor(i)=knownPa/rmsValueCal; %multiply measured values with this in order to get calibrated Pa

   
   %    
%    % End calibration
%    calfile = Dmeta(i).CalibrationFileEnd;
%    dat_temp = detrend(audioread(calfile));
%    
%    if strcmp(Dmeta(i).CalibrationStopIndEnd,'end')
%        dat = dat_temp(Dmeta(i).CalibrationStartIndEnd:end);
%    else
%        dat = dat_temp(Dmeta(i).CalibrationStartIndEnd:Dmeta(i).CalibrationStopIndEnd);       
%    end
%    dat_inf = audioinfo(calfile);
%    caldata(i,2) = caltarget/rms(dat);
%    disp(['This value shoud be 26.6. Check: ',num2str(rms(dat)*caldata(i,2))])
end

%% Sound exposures
SpawnSeisHydrophoneDataTreatment_testRMSogMaxogFrek(Dmeta,Tmeta,calibrationfactor)
%SpawnSeisHydrophoneDataTreatment_testRMSogMaxogFrekNoise(Dmeta,Tmeta,calibrationfactor)
%SpawnSeisHydrophoneDataTreatment(Dmeta,Tmeta,calibrationfactor)
%SpawnSeisHydrophoneDataTreatment_testFilter(Dmeta,Tmeta,calibrationfactor)
% %% Make map of hydrophone placements
% figure
% m_proj('albers equal-area','long',[5 5+8/60],'lat',[60+5/60 60+8/60]);
% m_grid('box','fancy','tickdir','in');   
% m_gshhs_f('patch',[.7 .9 .7]);
% hold on
% for i=1:length(Dmeta)
%     m_plot(Dmeta(i).LONdeg + Dmeta(i).LONmin/60,Dmeta(i).LATdeg + Dmeta(i).LATmin/60,'*')
% end
% colormap(flipud(copper));



