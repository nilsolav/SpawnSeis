% Hydrophone measurements for the SpawnSeis exposure project

%% Metadata for each hydrophone 
[~,~,Hmeta_raw] = xlsread('SpawnSeisHydrophoneMetaData.xlsx');
Hmeta=cell2struct(Hmeta_raw(2:end,:),Hmeta_raw(1,:),2);

% Metadata for each hydrophone deployment
[~,~,Dmeta_raw] = xlsread('SpawnSeisHydrophoneDeploymentMetaData.xlsx');
Dmeta=cell2struct(Dmeta_raw(2:end,:),Dmeta_raw(1,:),2);

% Get metadata for the treatments
f = 'D:\DATA\S2021xxx_PH.U.SVERDRUP II[1007]\CRUISE_LOG\ACTIVITY\treatements.xlsx';
[~,~,Tmeta_raw] = xlsread(f);
Tmeta=cell2struct(Tmeta_raw(2:end,:),Tmeta_raw(1,:),2);

%% Get calibration data per deployment
for i=2%1:length(Dmeta)
   caltarget = 10^(Dmeta(i).CalibrationLevel/20-6); % Pa
   % Beginning calibration
   calfile = Dmeta(i).CalibrationFileBeginning;
   dat_temp = detrend(audioread(calfile));
   
   if strcmp(Dmeta(i).CalibrationStopIndBeginning,'end')
       dat = dat_temp(Dmeta(i).CalibrationStartIndBeginning:end);
   else
       dat = dat_temp(Dmeta(i).CalibrationStartIndBeginning:Dmeta(i).CalibrationStopIndBeginning);       
   end
   
   dat_inf = audioinfo(calfile);
   caldata(i,1) = caltarget/rms(dat);
   disp(['This value shoud be 25.1. Check: ',num2str(rms(dat)*caldata(i,1))])
   
   % End calibration
   calfile = Dmeta(i).CalibrationFileEnd;
   dat_temp = detrend(audioread(calfile));
   
   if strcmp(Dmeta(i).CalibrationStopIndEnd,'end')
       dat = dat_temp(Dmeta(i).CalibrationStartIndEnd:end);
   else
       dat = dat_temp(Dmeta(i).CalibrationStartIndEnd:Dmeta(i).CalibrationStopIndEnd);       
   end
   dat_inf = audioinfo(calfile);
   caldata(i,2) = caltarget/rms(dat);
   disp(['This value shoud be 25.1. Check: ',num2str(rms(dat)*caldata(i,2))])
end

%% Sound exposures
SpawnSeisHydrophoneDataTreatment(Dmeta,Tmeta)


%% Make map of hydrophone placements
figure
m_proj('albers equal-area','long',[5 5+8/60],'lat',[60+5/60 60+8/60]);
m_grid('box','fancy','tickdir','in');   
m_gshhs_f('patch',[.7 .9 .7]);
hold on
for i=1:length(Dmeta)
    m_plot(Dmeta(i).LONdeg + Dmeta(i).LONmin/60,Dmeta(i).LATdeg + Dmeta(i).LATmin/60,'*')
end
colormap(flipud(copper));


