function SpawnSeisFileList(Tmeta_i,calibrationfactor)
% Extract time-pressure data for all sensors from one treatment
%

% Loop over the different deplyments relevant for this treatment
for j=1%:length(Tmeta_i.Hydrophone)
    % List all files from this deplyment
    fs=filesep;
    filelist = dir([Tmeta_i.Hydrophone(j).DataDir,fs,'*wav']);
    % Loop over files and get the start time vector for each file in the
    % deployment
    for k = 1:length(filelist)
        D = regexp(filelist(k).name,'\d*','Match');
        tmp2 = [D{2},'T',D{3}];
        T(k)=datenum(tmp2,'yyyymmddTHHMMSS');
    end
    % Get the file list for the treatment
    T0 = datenum(Tmeta_i.Start,'dd.mm.yyyy HH:MM:SS');
    T1 = datenum(Tmeta_i.End,'dd.mm.yyyy HH:MM:SS');
    ind = find((T > T0) & (T < T1)); % misses parts of the first file.
    if ind(1)>1
        ind = (min(ind)-1):1:max(ind); % Add the first file
    else
        ind = (min(ind)):1:max(ind); % But not if it already is the first file
    end
    % Extract the data for the treatment
    Dat.Pressure=[];
    Dat.Time=[];
    Dat.T0 = T0;
    % Loop over data files for this treatment
    for k=1:length(ind)
        F = fullfile(filelist(k).folder,filelist(k).name);
        ainf = audioinfo(F);
        dt = 1/ainf.SampleRate;
        dat_sub = audioread(F,'native');
        warning('Add calibration!')
        tim_sub = ((1:length(dat_sub))-1)*dt;% seconds
        % Time for this file in seconds relative to start of treatment (T0)
        T0_sub =  tim_sub + (T(ind(k))-T0)*3600;
        % This needs refactoring:
        Dat.Pressure = [Dat.Pressure dat_sub'];
        Dat.Time = [Dat.Time T0_sub];
    end
    tmpfil = ['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.mat'];
    save(tmpfil,'Dat','-v7.3')
    clear Dat
end
end