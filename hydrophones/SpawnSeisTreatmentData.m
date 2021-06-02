function SpawnSeisTreatmentData(Tmeta_i,calibrationfactor,tempdir)
% Extract time-pressure data for all sensors from one treatment
%

% Loop over the different deplyments relevant for this treatment
for j=1:length(Tmeta_i.Hydrophone)
    tmpfil = fullfile(tempdir,['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.mat']);
    if ~exist(tmpfil)
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
        if ~isempty(ind) % Check if there are missing times for the deployment
            if ind(1)>1
                ind = (min(ind)-1):1:max(ind); % Add the first file
            else
                ind = (min(ind)):1:max(ind); % But not if it already is the first file
            end
            % Get calibration factor
            cal = calibrationfactor(Tmeta_i.Hydrophone(j).Dmeta_index,1);
            
            % Extract the data for the treatment
            Dat.Pressure=zeros(1,400000000);
            Dat.Time=zeros(1,400000000);
            Dat.T0 = T0;
            k1=1;
            % Loop over data files for this treatment
            for k=1:length(ind)
                try
                    F = fullfile(filelist(ind(k)).folder,filelist(ind(k)).name);
                    ainf = audioinfo(F);
                    dt = 1/ainf.SampleRate;
                    dat_sub = audioread(F,'native');
                    tim_sub = ((1:length(dat_sub))-1)*dt;% seconds
                    % Time for this file in seconds relative to start of treatment (T0)
                    T0_sub =  tim_sub + (T(ind(k))-T0)*3600*24;
                    % Adding data
                    k2 = k1-1+length(dat_sub);
                    Dat.Pressure(k1:k2) = double(dat_sub')*cal;
                    Dat.Time(k1:k2) = T0_sub;
                    k1 = k2+1;
                catch
                    warning(['Cannot read ',F,' -> skipping file'])
                end
            end
            % Trim data files
            Dat.Pressure = Dat.Pressure(1:k2);
            Dat.Time = Dat.Time(1:k2);
            
            save(tmpfil,'Dat','-v7.3')
            clear Dat
        else
            warning(['No data for Block',num2str(Tmeta_i.BlockNo),' Treatment',num2str(Tmeta_i.TreatmentNo),' Hydrophone',num2str(j)])
        end
    end
end
end