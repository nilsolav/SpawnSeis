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
        
        % Detect and filter pulses (Nils Olav)
        if ~exist(pulsefil)
            p = Dat.Pressure;
            t = Dat.Time;
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
hold on
for i=1:length(loc)
    indpulse = t > (loc(i)-par.tmin) & t < (loc(i)+par.tmax);
    % Adding ten samples to be on the safe side. Sloppy, but works.
    if (length(t(indpulse))+10)>(par.tmin+par.tmax)*par.Fs
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

function pulse=AnalyzePulse(t,p,t0,plt,frek)
% Analyze each individual pulse
figure
plot(t,p)

%sampling rate
dt=t(3)-t(2);
Fs=round(1/dt);

%find index for peak at t0
ind_tid=find(t>=t0,1)

%select one second around peak
tid=ind_tid-round(0.3*Fs):ind_tid+round(0.7*Fs);

%Do the same ting for noise
tidN=tid-round(2*Fs);

%Select 1 second long signal to analyze
S1=p(tid);%signal
S1N=p(tidN);%st�y

%Find and save max and min:
[p1_max,t1_max]=max(S1);
[p1_min,t1_min]=min(S1);

[p1N_max,t1N_max]=max(S1N);
[p1N_min,t1N_min]=min(S1N);

pulse.pospeakpressure=p1_max;
pulse.negpeakpressure=p1_min;

pulse.pospeakpressureN=p1N_max;
pulse.negpeakpressureN=p1N_min;

%Find and save Ex and SEL:

pulse.Ex=dt*sum(S1.^2); %energi i signal
pulse.ExN=dt*sum(S1N.^2); %same som rms
pulse.SEL=10*log10(pulse.Ex/1e-12)
pulse.SELN=10*log10(noise.Ex/1e-12)

%% frekvensanalyse
if frek==1;
    tuk=tukeywin(length(S1),0.3); %Tapering: lagar vindu som gir ein glatt overgang ved � setje start og sluttverdi p� tidsvindu til 0
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
    
    pulse.ESD=ESD;
    pulse.ESDN=ESDN;
    pulse.F=frekv
    
 %% plot   
    if plt==1
       figure
    subplot(2,1,1)
    plot(t,p)
    hold on
    plot(t(tid),S1,'k')
    plot(t(tidN),S1N,'r')
   % legend(['max peak=' num2str(round(20*log10(pulse.pospeakpressure/1e-6),1)) ' dB re 1 \muPa'],['1 sek seismic SEL=' num2str(round(pulse.SEL),1)) 'dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round((pulse.SELN),1) 'dB re 1 \muPa^2\cdots'])
    %title(['1 second ' tittel], 'Interpreter', 'none')
    xlabel('clock time, UTC')
    ylabel('Pa')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
    set(gca,'fontsize', 22);
    
    subplot(2,1,2)

    plot(frek,10*log10(ESD/1e-12),'k')
    hold on
    plot(frek,10*log10(ESDN/1e-12),'r')
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
    legend(['max peak=' num2str(round(20*log10(p0/1e-6),1)) ' dB re 1 \muPa'],['1 sek seismic SEL=' num2str(round(10*log10(SEL_t/1e-12),1)) 'dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(10*log10(SEL_tN/1e-12),1)) 'dB re 1 \muPa^2\cdots'])
    %title(['1 second ' tittel], 'Interpreter', 'none')
    xlabel('clock time, UTC')
    ylabel('Pa')
    set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
    set(gca,'fontsize', 22);
    
    subplot(2,1,2)
    plot(frek,P1,'k')
    hold on;
    plot(frek,P1N,'r')
    legend('seismic','ambient')
    ylabel('Pa')
    xlabel('Frequency, Hz')
    
    %  plot(frek,10*log10(ESD/1e-12),'k')
    %  hold on
    %   plot(frek,10*log10(ESDN/1e-12),'r')
    %      title('Energy Spectral Density')
    % ylabel({'dB re 1 \muPa^2\cdots/Hz'})
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