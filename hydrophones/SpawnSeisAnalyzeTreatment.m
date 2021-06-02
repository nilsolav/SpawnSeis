function Pulses = SpawnSeisAnalyzeTreatment(Tmeta_i,tempdir,par)


% Extract time-pressure data for all sensors from one treatment
%

% Loop over the different deplyments relevant for this treatment
for j=1:length(Tmeta_i.Hydrophone)
    % read TEMP files
    tmpfil = ['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.mat'];
    figfil = ['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.png'];
    load(fullfile(tempdir,tmpfil),'Dat'); % Loads DAT
    % Dat.Pressure
    % Dat.Time
    % Dat.T0
    
    % Filter data prior to analysis
    %par.tmin = 4;
    %par.tmax = 4;
    
    % Detect and filter pulses (Nils Olav)
    Pulses = DetectPulse(Dat);
    
    % Analyze individual pulses 
    for k=1:length(Pulses)
        % Select pulse
        ind = Dat.Time>(Pulses(k).t0-par.tmin) & Dat.Time<(Pulses(k).t0+par.tmax);
        t = Dat.Time(ind);
        p = Dat.Pressure(ind);
        
        % Calculate stats per pulse (Tonje)
        Pulses(k).pulse = AnalyzePulse(t,p,Pulses(k).t0);
    end
    
    
%     % Plot files in headless mode
%     f = figure('visible', 'off');
%     plot(Dat.Time,upper,'k',Dat.Time,lower,'k')
%     title(Tmeta_i.Comment)
%     xlabel('Time relative to start Treatment (s)')
%     ylabel('Pressure envelope (Pa)')
%     print(fullfile(tempdir,figfil),'-dpng')
%     close(f)

end
end

function Pulses = DetectPulse(Dat)
% detects complete pulses (time and peakpressures)
Pulses(1).t0=95.36;
Pulses(1).pp=43;

end

function pulse=AnalyzePulse(t,p,t0)
% Analyze each individual pulse
figure
plot(t,p)
pulse=NaN;
end