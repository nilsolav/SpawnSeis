function Pulses = SpawnSeisAnalyzeTreatment(Tmeta_i,tempdir,par)


% Extract time-pressure data for all sensors from one treatment
%

% Loop over the different deplyments relevant for this treatment
for j=1:length(Tmeta_i.Hydrophone)
    % read TEMP files
    tmpfil = fullfile(tempdir,['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.mat']);
    figfil = fullfile(tempdir,['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.png']);
    pulsefil = fullfile(tempdir,['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'pulses.mat']);
    if exist(tmpfil)
        load(tmpfil,'Dat'); % Loads DAT
        
        % Filter pressure data
        [C,D] = butter(3,[5 10000]/(par.Fs/2),'bandpass');
        p = filtfilt(C,D,Dat.Pressure);% ca 30 sec run time
        t = Dat.Time;
        
        % Detect pulses
        if ~exist(pulsefil)
            DetectedPulses = DetectPulse(t,p,par);
            save(pulsefil,'DetectedPulses')
        else
            load(pulsefil,'DetectedPulses')
        end
        
        % Plot files in headless mode
        f = figure('visible', 'off');
        plot(DetectedPulses.t0_f,DetectedPulses.pp_f,'r.',DetectedPulses.t0,DetectedPulses.pp,'b.')
        title(Tmeta_i.Hydrophone(j).Comment)
        xlabel('Time relative to start Treatment (s)')
        ylabel('Positive Peak Pressure (Pa)')
        legend('Incomplete pulses','Valid pulses')
        print(figfil,'-dpng')
        close(f)
        
        
        % Analyze individual pulses
        %     for k=1:length(DetectedPulses.pp)
        %         % Select pulse
        %         ind = Dat.Time>(DetectedPulses.t0(k)-par.tmin) & Dat.Time<(DetectedPulses.t0(k)+par.tmax);
        %         t = Dat.Time(ind);
        %         p = Dat.Pressure(ind);
        %
        %         % Calculate stats per pulse (Tonje)
        %         frek=1;%include frequency analysis
        %         plt=1; %plot figures for pulse analysis
        %         Pulses(k).pulse = AnalyzePulse(t,p,Pulses(k).t0,plt,frek);
        %     end
    end
end
end

function Pulses = DetectPulse(t,p,par)
% Detects the pulses (time and positive peak pressures)

% Find peak candidates
[pks,loc]=findpeaks(p,t,'MinPeakDistance',par.minpeakdistance);

% Loop and remove the peaks closer than par.tmin/par.tmax to the break
locind = false(1,length(loc));
% Add variable for the negative pulses
pospeakpressure = NaN(size(loc));
negpeakpressure = NaN(size(loc));
pospeakpressureN = NaN(size(loc));
negpeakpressureN = NaN(size(loc));
Ex = NaN(size(loc));
ExN = NaN(size(loc));
SEL = NaN(size(loc));
SELN = NaN(size(loc));

hold on
for i=1:length(loc)
    % Get minimum pressure
    indpulse = t > (loc(i)+par.tmin) & t < (loc(i)+par.tmax);
    % Adding ten samples to be on the safe side. Sloppy, but works.
    if (length(t(indpulse))+10)>(-par.tmin+par.tmax)*par.Fs
       locind(i)= true;
    else
        locind(i)=false;
    end
    
    % Analyze each pulse
    dum = AnalyzePulse(t(indpulse),p(indpulse),loc(i),par,false,false);
    pospeakpressure(i) = dum.pospeakpressure;
    negpeakpressure(i) = dum.negpeakpressure;
    pospeakpressureN(i) = dum.pospeakpressureN;
    negpeakpressureN(i) = dum.negpeakpressureN;
    Ex(i) = dum.Ex;
    ExN(i) = dum.ExN;
    SEL(i) = dum.SEL;
    SELN(i) = dum.SELN;
end

% testplot
%plot(loc(locind),pks(locind),'b-*',loc(~locind),pks(~locind),'r*')
%semilogy(loc(locind),pks(locind),'b-*',loc(~locind),pks(~locind),'r*')
%semilogy(loc(locind),pks(locind),'b*')%,loc(~locind),pks(~locind),'r*')
%plot(t,p,'k',loc(locind),pks(locind),'b*',loc(~locind),pks(~locind),'r*')

Pulses.t0 = loc;
Pulses.ind_ok = locind;
Pulses.pospeakpressure = pospeakpressure;
Pulses.negpeakpressure = negpeakpressure;
Pulses.pospeakpressureN = pospeakpressureN;
Pulses.negpeakpressureN = negpeakpressureN;
Pulses.Ex = Ex;
Pulses.ExN =ExN;
Pulses.SEL = SEL;
Pulses.SELN = SELN;


end

function pulse=AnalyzePulse(t,p,t0,par,plt,frek)
% Analyze each individual pulse

%sampling rate
Fs = par.Fs;
dt=1/Fs;

%find index for peak at t0
ind_tid=find(t>=t0,1);

% Time relative to peak
t_rel = t-t(ind_tid);

%Select one second around peak
tid = t_rel>par.SELinterval(1) & t_rel<par.SELinterval(2);

%Do the same ting for noise before peak
tidN = t_rel>par.noiseStart(1) & t_rel<par.noiseStart(2);

%Select 1 second long signal to analyze
S1=p(tid);%signal
S1N=p(tidN);%støy

%Find and save max and min:
[p1_max,~]=max(S1);
[p1_min,~]=min(S1);

[p1N_max,~]=max(S1N);
[p1N_min,~]=min(S1N);

% Initialize
pulse.pospeakpressure=NaN;
pulse.negpeakpressure=NaN;
pulse.pospeakpressureN=NaN;
pulse.negpeakpressureN=NaN;
pulse.Ex=NaN;
pulse.ExN=NaN; 
pulse.SEL=NaN;
pulse.SELN=NaN;

% Add values
pulse.pospeakpressure=p1_max;
pulse.negpeakpressure=p1_min;

pulse.pospeakpressureN=p1N_max;
pulse.negpeakpressureN=p1N_min;

%Find and save Ex and SEL:

pulse.Ex=dt*sum(S1.^2); %energi i signal
pulse.ExN=dt*sum(S1N.^2); %same som rms
pulse.SEL=10*log10(pulse.Ex/1e-12);
pulse.SELN=10*log10(pulse.ExN/1e-12);

%% frequency analysis
if frek
    tuk=tukeywin(length(S1),0.3)'; %Tapering: lagar vindu som gir ein glatt overgang ved å setje start og sluttverdi på tidsvindu til 0
    S=tuk.*S1; %1 sekund signal-sekvens
    %tek fft;
    SN=tuk.*S1N;
    Y=fft(S); %1 sekund signal-sekvens
    L=length(tid(1,:));
    P2=abs(Y/L);
    P1=P2(1:L/2+1); %tek halve spekteret
    P1(2:end-1)=2*P1(2:end-1); %gongar mesteparten av spekteret med 2
    frekv=Fs*(0:(L/2))/L;
    ESD=((abs(P1)).^2)*(L/(2*Fs));
    
    YN=fft(SN);
    P2N=abs(YN/L);
    P1N=P2N(1:L/2+1);
    P1N(2:end-1)=2*P1N(2:end-1);
    ESDN=((abs(P1N)).^2)*(L/(2*Fs));
    
    pulse.ESD=10*log10(ESD/1e-12);
    pulse.ESDN=10*log10(ESDN/1e-12);
    pulse.F=frekv;
    
    %% plot
    if plt
        figure
        subplot(2,1,1)
        plot(t,p)
        hold on
        plot(t(tid),S1,'k')
        plot(t(tidN),S1N,'r')
        legend(['max peak=' num2str(round(20*log10(pulse.pospeakpressure/1e-6),1)) ' dB re 1 \muPa'],['1 sek seismic SEL=' num2str(round(pulse.SEL,1)) ' dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(pulse.SELN,1)) ' dB re 1 \muPa^2\cdots'])
        %title(['1 second ' tittel], 'Interpreter', 'none')
        xlabel('clock time, UTC')
        ylabel('Pa')
        set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
        set(gca,'fontsize', 22);
        
        subplot(2,1,2)
        
        plot(pulse.F,pulse.ESD,'k')
        hold on
        plot(pulse.F,pulse.ESDN,'r')
        title('Energy Spectral Density')
        ylabel({'dB re 1 \muPa^2\cdots/Hz'})
        legend('seismic','ambient')
        set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
        xlim([0 1000])
        
        xlabel('Frequency, Hz')
        set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
            'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
        set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
        
        
        figure
        subplot(2,1,1)
        plot(t,p)
        hold on
        plot(t(tid),S,'k')
        plot(t(tidN),SN,'r')
        legend(['max peak=' num2str(round(20*log10(pulse.pospeakpressure/1e-6),1)) ' dB re 1 \muPa'],['1 sek seismic SEL=' num2str(round(pulse.SEL,1)) ' dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(pulse.SELN,1)) ' dB re 1 \muPa^2\cdots'])
        
        %title(['1 second ' tittel], 'Interpreter', 'none')
        xlabel('clock time, UTC')
        ylabel('Pa')
        set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
        set(gca,'fontsize', 22);
        
        subplot(2,1,2)
        plot(pulse.F,P1,'k')
        hold on;
        plot(pulse.F,P1N,'r')
        legend('seismic','ambient')
        ylabel('Pa')
        xlabel('Frequency, Hz')
        
        legend('seismic','ambient')
        set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
        xlim([0 1000])
        set(gca,'fontsize', 22);
        
        set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
            'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
        set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
        
        
        
    end
    
end

end      