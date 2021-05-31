% Hydrophone measurements for the SpawnSeis exposure project
% Requires the m_map package

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

% Metadata for each hydrophone 
[~,~,Hmeta_raw] = xlsread('SpawnSeisHydrophoneMetaData.csv');
Hmeta=cell2struct(Hmeta_raw(2:end,:),Hmeta_raw(1,:),2);

% Metadata for each hydrophone deployment
[~,~,Dmeta_raw] = xlsread('SpawnSeisHydrophoneDeploymentMetaData.csv');
Dmeta=cell2struct(Dmeta_raw(2:end,:),Dmeta_raw(1,:),2);

% Get metadata for the treatments
[~,~,Tmeta_raw] = xlsread('treatments.csv');
Tmeta=cell2struct(Tmeta_raw(2:end,:),Tmeta_raw(1,:),2);

%% Check if calibration files exists
% excluding sound trap and vessel mounted hydrophones, i.e. i \in [3 4 8 9 10 11 12 13])
for i = [1 2 5 6 7 14 15 17 18 19 20]
    calfile_e = fullfile(rootdir,Dmeta(i).StartTime(7:10),Dmeta(i).CalibrationFileEnd);
    calfile_b = fullfile(rootdir,Dmeta(i).StartTime(7:10),Dmeta(i).CalibrationFileBeginning);
    if ~exist(calfile_b)
        disp(['Index:',num2str(i),'; Missing: ',calfile_b])
    end
    if ~exist(calfile_e)
        disp(['Index:',num2str(i),'; Missing: ',calfile_e])
    end
end

%% Get calibration data per deployment and store in .mat file
for i = [1 2 5 6 7 14 15 17 18 19 20]
   knownPa = 10^((Dmeta(i).CalibrationLevel/20)-6); % Pa reference pressure rms-value in calibrator with naxys coupler (Se doc in H2... folder)
   % Beginning calibration
   Calfile{1} = fullfile(rootdir,Dmeta(i).StartTime(7:10),Dmeta(i).CalibrationFileBeginning);
   % End calibration
   Calfile{2} = fullfile(rootdir,Dmeta(i).StartTime(7:10),Dmeta(i).CalibrationFileEnd);
   % Loop over beginning and end calibration
   for j=1:2
       disp([i j])
       calfile = Calfile{j};
       % Test if calibratiuon file exist
       if exist(calfile)
           % Read calibration files and convert to Pa
           dat_temp = audioread(calfile,'native'); %les inn rå-versjonen, utan native blir normalisert til 1.
           dat = detrend(double(dat_temp));
           dat_inf = audioinfo(calfile);
           caldata(i,1) = knownPa/rms(dat);
           disp(['This value should be 26.6. Check: ',num2str(rms(dat)*caldata(i,1))])
           Fs = dat_inf.SampleRate;
           rmsValueCal=rms(dat);
           calibrationfactor(i,j)=knownPa/rmsValueCal; %multiply measured values with this in order to get calibrated Pa
       else
            disp('This value should be 26.6. Check: NaN')
           calibrationfactor(i,j) = NaN;
       end
   end
end
save calibrationfactor.mat calibrationfactor

%% Append Dmeta with file placement for raw data
for i=1:length(Tmeta)
%    disp([datestr(datenum(Tmeta(i).Start),'dd.mm.yyyy HH:MM:SS'),'   ',  datestr(datenum(Tmeta(i).End,'dd.mm.yyyy HH:MM:SS'))])
%    disp(' ')
    Tmeta(i).DataDir = NaN;
    % Relevant deplyments for this treatment
    k=1;
    for j=[1 2 5 6 7 14 15 17 18 19 20]
%        disp([datestr(datenum(Dmeta(j).StartTime,'dd.mm.yyyy HH:MM:SS')),'   ',  datestr(datenum(Dmeta(j).StopTime,'dd.mm.yyyy HH:MM:SS'))])
        if datenum(Tmeta(i).Start,'dd.mm.yyyy HH:MM:SS') > datenum(Dmeta(j).StartTime,'dd.mm.yyyy HH:MM:SS') && datenum(Tmeta(i).End,'dd.mm.yyyy HH:MM:SS') < datenum(Dmeta(j).StopTime,'dd.mm.yyyy HH:MM:SS')
            Tmeta(i).Hydrophone(k).DataDir = fullfile(rootdir,Dmeta(j).StartTime(7:10),[Dmeta(j).CalibrationFileBeginning(1:8),'_PH.U.SVERDRUP II[1007]'],'EXPERIMENTS','HYDROPHONES',Dmeta(j).Folder);
            Tmeta(i).Hydrophone(k).Dmeta_index = j;
            Tmeta(i).Hydrophone(k).Comment = Dmeta(j).Comment;
            k=k+1;
        end        
    end    
end

%% Analyze each treatment
SpawnSeisHydrophoneDataTreatment_all(Dmeta,Tmeta,calibrationfactor)

%% Sound exposures
SpawnSeisHydrophoneDataTreatment_testRMSogMaxogFrek(Dmeta,Tmeta,calibrationfactor)
%SpawnSeisHydrophoneDataTreatment_testRMSogMaxogFrekNoise(Dmeta,Tmeta,calibrationfactor)
%SpawnSeisHydrophoneDataTreatment(Dmeta,Tmeta,calibrationfactor)
%SpawnSeisHydrophoneDataTreatment_testFilter(Dmeta,Tmeta,calibrationfactor)
%% Make map of hydrophone placements
figure
m_proj('albers equal-area','long',[5 5+8/60],'lat',[60+5/60 60+8/60]);
m_grid('box','fancy','tickdir','in');   
m_gshhs_f('patch',[.7 .9 .7]);
hold on
for i=1:length(Dmeta)
    m_plot(Dmeta(i).LONdeg + Dmeta(i).LONmin/60,Dmeta(i).LATdeg + Dmeta(i).LATmin/60,'*')
%    m_text(Dmeta(i).LONdeg + Dmeta(i).LONmin/60,Dmeta(i).LATdeg + Dmeta(i).LATmin/60,'*')
end
colormap(flipud(copper));



