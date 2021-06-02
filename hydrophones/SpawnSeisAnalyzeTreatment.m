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
    DetectedPulses = DetectPulse(Dat);
    
    % Analyze individual pulses 
    for k=1:length(DetectedPulses.pp)
        % Select pulse
        ind = Dat.Time>(DetectedPulses.t0(k)-par.tmin) & Dat.Time<(DetectedPulses.t0(k)+par.tmax);
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
% Detects the pulses (time and positive peak pressures)

% Subset for testing
% Pulses(1).t0=95.36;
% Pulses(1).pp=43;
% par.tmin = 30;
% par.tmax = 30;
% ind = Dat.Time>(Pulses(1).t0-par.tmin) & Dat.Time<(Pulses(1).t0+par.tmax);
% p = Dat.Pressure(ind);
% t = Dat.Time(ind);
%  ind = Dat.Time>6670 & Dat.Time<6720;
%  p = Dat.Pressure(ind);
%  t = Dat.Time(ind);

p = Dat.Pressure;
t = Dat.Time;

% par.tmin = 4;
% par.tmax = 4;
% Fs = par.Fs;

plot(t,p)

% Find peak candidates
[pks,loc]=findpeaks(p,t,'MinPeakDistance',par.minpeakdistance);

% Loop and remove the peaks closer than par.tmin/par.tmax to the break
locind = false(1,length(loc));
hold on
for i=1:length(loc)
    indpulse = t > (loc(i)-par.tmin) & t < (loc(i)+par.tmax);
    % Adding ten samples to be on the safe side. Sloppy, but works.
    if (length(t(indpulse))+10)>par.minpeakdistance*Fs
       locind(i)= true;
    else
        locind(i)=false;
    end
end
% testplot
%plot(loc(locind),pks(locind),'b-*',loc(~locind),pks(~locind),'r*')
%semilogy(loc(locind),pks(locind),'b-*',loc(~locind),pks(~locind),'r*')
%semilogy(loc(locind),pks(locind),'b*')%,loc(~locind),pks(~locind),'r*')
%plot(t,p,'k',loc(locind),pks(locind),'b*',loc(~locind),pks(~locind),'r*')

Pulses.t0 = loc(locind);
Pulses.pp = pks(locind);
Pulses.t0_f = loc(~locind);
Pulses.pp_f = pks(~locind);

end

function pulse=AnalyzePulse(t,p,t0)
% Analyze each individual pulse
figure
plot(t,p)
pulse=NaN;
end