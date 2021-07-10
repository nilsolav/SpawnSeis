function Pulses = SpawnSeisAnalyzeTreatment(Tmeta_i,tempdir,par)


% Extract time-pressure data for all sensors from one treatment

% Loop over the different deplyments relevant for this treatment
for j=1:length(Tmeta_i.Hydrophone)
    % read TEMP files
    tmpfil = fullfile(tempdir,['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_Hydr',num2str(j),'.mat']);
    
    % Generate descriptive figure file names
    figfil = fullfile(tempdir,['Block',num2str(Tmeta_i.BlockNo),'_Treat',num2str(Tmeta_i.TreatmentNo),'_',Tmeta_i.Treatment,...
        '_',Tmeta_i.Hydrophone(j).DeplNumber,'_Location_',Tmeta_i.Hydrophone(j).Location]);

    % Temporary pulse file to avoid detecting the pulses at each run
    pulsefil = [figfil,'_pulses.mat'];
    
    % Generate a separate folder for the pulses
    if ~exist(figfil)
        mkdir(figfil)
    end
    
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
        D = DetectedPulses;
        f = figure('visible', 'off');
        plot(D.t0(D.ind_ok)/60,D.pospeakpressure(D.ind_ok),'k.',...
            D.t0(D.ind_ok)/60,D.negpeakpressure(D.ind_ok),'r.')
        title([Tmeta_i.Treatment ', ' Tmeta_i.Hydrophone(j).Comment])
        legend({'Positive peak','Negative Peak'})
        xlabel('Time relative to start Treatment (min)')
        ylabel('Peak Pressure (Pa)')
        print([figfil,'_peakpressure'],'-dpng')
        close(f)
        
        f = figure('visible', 'off');
        plot(D.t0(D.ind_ok)/60,D.SEL(D.ind_ok),'k.',...
            D.t0(D.ind_ok)/60,D.SELN(D.ind_ok),'r.')
        title([Tmeta_i.Treatment ', ' Tmeta_i.Hydrophone(j).Comment])
        legend({'Pulse','Background Noise'})
        xlabel('Time relative to start Treatment (min)')
        ylabel('Sound Exposure Level (SEL) (dB re 1\muPa^2 s)')
        print([figfil,'_SEL'],'-dpng')
        close(f)
        
        % Estimate SEL average
        D.Ex_all=(1/par.Fs)*sum(p.^2);
        D.SEL_all_dB=10*log10(D.Ex_all/1e-12);
        D.Ex_all_1sec=D.Ex_all/((t(end)-t(1))*22/30); % (tilnærming: lydopptak 22/30 delar av tida)
        D.SEL_all_1sec_dB=10*log10(D.Ex_all_1sec/1e-12);
        D.rms=20*log10(rms(p)./1e-6);
        
        f = figure('visible', 'off');
        plot(t/60,p,'k')
        legend(['cumSEL=' num2str(round(D.SEL_all_dB,1)) 'dB re 1\muPa^2 s,' 10 'rms=' num2str(round(D.rms,1)) ' dB re 1\muPa'])
        title([Tmeta_i.Treatment ', ' Tmeta_i.Hydrophone(j).Comment])
        xlabel('Time relative to start Treatment (min)')
        ylabel('Pa')
        print([figfil,'_raw'],'-dpng')
        close(f)
        
        % Write data from plots to csv file
        Header = {'t0','pospeakpressure','negpeakpressure','SEL','SELN'};
        writecell(Header,[figfil,'.csv'],'Delimiter',';')
        Data = [D.t0(D.ind_ok)/60;D.pospeakpressure(D.ind_ok);D.negpeakpressure(D.ind_ok);D.SEL(D.ind_ok);D.SELN(D.ind_ok)]';
        writematrix(Data,[figfil,'.csv'],'Delimiter',';','WriteMode','append')
        
        % Pick out pulses for more detailed analyses
        peak_ok=D.pospeakpressure(D.ind_ok);
        t_ok=D.t0(D.ind_ok);
        [~,indp(1)]=max(peak_ok);
        [~,indp(2)]=min(peak_ok);
        
        % In addition to max and mean some more peaks are selected based 
        % on the percentile
        tempd= cumsum(peak_ok);
        tempd=tempd./max(tempd);  % cumulative sum normalized to 1
        indp(3)=max(find(tempd<par.prctile(1)));
        indp(4)=max(find(tempd<par.prctile(2)));
        indp(5)=max(find(tempd<par.prctile(3)));
        
        % Analyse each selected pulse
        for b=1:5
            indpulse = t > (t_ok(indp(b))+par.tmin) & (t < t_ok(indp(b))+par.tmax);
            dum = AnalyzePulse(t(indpulse),p(indpulse),t_ok(indp(b)),par,true,true,figfil);
        end
    end
end
end

function Pulses = DetectPulse(t,p,par)
% Detects the pulses (time and positive peak pressures)

% Find peak candidates
[~,loc]=findpeaks(p,t,'MinPeakDistance',par.minpeakdistance);

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
textprogressbar('Analyze pulses: ')
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
    textprogressbar(100*(i/length(loc)))
end
textprogressbar(' Pulse analysis completed.')

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


function textprogressbar(c)
% This function creates a text progress bar. It should be called with a 
% STRING argument to initialize and terminate. Otherwise the number correspoding 
% to progress in % should be supplied.
% INPUTS:   C   Either: Text string to initialize or terminate 
%                       Percentage number to show progress 
% OUTPUTS:  N/A
% Example:  Please refer to demo_textprogressbar.m
% Author: Paul Proteus (e-mail: proteus.paul (at) yahoo (dot) com)
% Version: 1.0
% Changes tracker:  29.06.2010  - First version
% Inspired by: http://blogs.mathworks.com/loren/2007/08/01/monitoring-progress-of-a-calculation/
%% Initialization
persistent strCR;           %   Carriage return pesistent variable
% Vizualization parameters
strPercentageLength = 10;   %   Length of percentage string (must be >5)
strDotsMaximum      = 10;   %   The total number of dots in a progress bar
%% Main 
if isempty(strCR) && ~ischar(c),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(c),
    % Progress bar - initialization
    fprintf('%s',c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c),
    % Progress bar  - termination
    strCR = [];  
    fprintf([c '\n']);
elseif isnumeric(c)
    % Progress bar - normal progress
    c = floor(c);
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(c/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end
end
function pulse=AnalyzePulse(t,p,t0,par,plt,frek,figdir)


if nargin<7
    figdir= '.';
end
% Analyze each individual pulse

%sampling rate
Fs = par.Fs;
dt=1/Fs;

%find index for peak at t0
ind_tid=find(t>=t0,1);

% Time relative to peak
t_rel = t-t(ind_tid);

%Select indexes for one second around peak
tid = t_rel>par.SELinterval(1) & t_rel<=par.SELinterval(2);

%Do the same ting for noise before peak
tidN = t_rel>par.noiseStart(1) & t_rel<=par.noiseStart(2);

%Select 1 second long signal to analyze
S1=p(tid);%signal  % Må bli 1 sekund = Fs antall samples
S1N=p(tidN);%støy

% if (length(S1)<Fs || length(S1)>Fs)
%      warning(['Not 1 second around pulse for analysis, length of S1: ' num2str(length(S1))])
%
% elseif length(S1N) <Fs || length(S1N) > Fs
%     warning(['Not 1 second around noise for analysis, length of S1N: ' num2str(length(S1N))])
%
% else
%
% end
%

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

% Calculate signal values
pulse.pospeakpressure=p1_max;
pulse.negpeakpressure=p1_min;
pulse.Ex=dt*sum(S1.^2); %energi i signal
pulse.SEL=10*log10(pulse.Ex/1e-12);

% Calculate noise values (if interval is present in the pulse)
if length(find(tidN))>1
    pulse.pospeakpressureN=p1N_max;
    pulse.negpeakpressureN=p1N_min;
    pulse.ExN=dt*sum(S1N.^2); %same som rms
    pulse.SELN=10*log10(pulse.ExN/1e-12);
end

%% frequency analysis
if frek
    if length(S1)==48000 && length(S1N)==48000
        tuk=tukeywin(length(S1),0.3)'; %Tapering: lagar vindu som gir ein glatt overgang ved å setje start og sluttverdi på tidsvindu til 0
        S=tuk.*S1; %1 sekund signal-sekvens
        %tek fft;
        SN=tuk.*S1N;
        Y=fft(S); %1 sekund signal-sekvens
        L=length(S);
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
    end
    
    %% plot
    if plt
        
        s = get(0, 'ScreenSize');
        f = figure('Position', [0 0 s(3) s(4)],'visible', 'off');
        subplot(2,1,1)
        plot(t,p)
        hold on
        plot(t(tid),S1,'k')
        plot(t(tidN),S1N,'r')
        
        legend(['max peak=' num2str(round(20*log10(pulse.pospeakpressure/1e-6),1)) ' dB re 1 \muPa'],['1 sek signal SEL=' num2str(round(pulse.SEL,1)) ' dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(pulse.SELN,1)) ' dB re 1 \muPa^2\cdots'],'Location','SouthWest')
        %title(['1 second ' tittel], 'Interpreter', 'none')
        xlabel('seconds')
        ylabel('Pa')
        set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
        set(gca,'fontsize', 25);
        
        if isfield(pulse,'F') % Check if data exist before plotting
            subplot(2,1,2)
            plot(pulse.F,pulse.ESD,'k')
            hold on
            plot(pulse.F,pulse.ESDN,'r')
            title('Energy Spectral Density')
            ylabel({'dB re 1 \muPa^2\cdots/Hz'})
            legend('signal','ambient')
            set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
            xlim([0 1000])
            
            xlabel('Frequency, Hz')
            set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
            set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
        end
        print(f,fullfile(figdir,['Onesecond_ESD_' num2str(round(t0))]),'-dpng')
        close(f)
        
        f = figure('Position', [0 0 s(3) s(4)],'visible', 'off');
        subplot(2,1,1)
        plot(t,p)
        hold on
        plot(t(tid),S,'k')
        plot(t(tidN),SN,'r')
        legend(['max peak=' num2str(round(20*log10(pulse.pospeakpressure/1e-6),1)) ' dB re 1 \muPa'],['1 sek signal SEL=' num2str(round(pulse.SEL,1)) ' dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(pulse.SELN,1)) ' dB re 1 \muPa^2\cdots'])
        
        %title(['1 second ' tittel], 'Interpreter', 'none')
        xlabel('seconds')
        ylabel('Pa')
        set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
        set(gca,'fontsize', 25);
        
        subplot(2,1,2)
        plot(pulse.F,P1,'k')
        hold on;
        plot(pulse.F,P1N,'r')
        
        ylabel('Pa')
        xlabel('Frequency, Hz')
        
        legend('signal','ambient')
        
        xlim([0 1000])
        set(gca,'fontsize', 25);
        
        set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',25, ...
            'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
        set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
        
        print(f,fullfile(figdir,['Onesecond_Pa_' num2str(round(t0))]),'-dpng')
        close(f)
        
    end
    
end

end