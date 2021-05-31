function plotSamples(p,t,t0)
% 
% Analyze each pulse as a function of time
% 
% Input:
% p - pressure
% t - time
% t0 - time at pulse
% 
% Output:
% pulse.pospeakpressure [Pa]
% pulse.negpeakpressure [Pa]
% pulse.Ex (sound exposure, T=1s, [\muPa^2 s])
% pulse.SEL (sound exposure level, T=1s, [dB re (1 \muPa)^2 (1 s)])
%
% pulse.ESD : energy spectral density
% pulse.F : frequency vector
% 

keyboard
Fs=1/dt;
MPH=max(data)*0.6;    

     MPD=0.9; %minimum peak distance, 9.95 seconds
     l=length(data);
   
        [PKS,LOCS] = findpeaks(abs(data),l,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
  [PKS,LOCS] = findpeaks(abs(data),'MinPeakHeight',MPH); %tar vekk sikkerhetssone
    if ~isempty(PKS)
       if length(PKS)>1
           [PKS, pos]=max(PKS);
           LOCS=LOCS(pos);
                      
       end
    end
    
                       
                        

                        
                        tid=round(LOCS)-15000:1:round(LOCS+Fs)-15000-1;%1 sekund signal
                        tidN=tid+8*Fs;
                        
                        if length(tidsakse)<tidN(end)
                          tidN=tid-2*Fs;  
                        end
                        % 1 sekund sekvens
                        H=data;
                        S1=H(tid);
                        S1N=H(tidN);
                        tuk=tukeywin(length(S1),0.3); %lagar vindu som gir ein glatt overgang ved å setje start og sluttverdi på tidsvindu til 0
                        S=tuk.*S1; %1 sekund signal-sekvens
                        %tek fft;
                        SN=tuk.*S1N;
                        
              
                        SEL_t=dt*sum(S1.^2); %energi i signal 
                        SEL_tN=dt*sum(S1N.^2);
                        
                        
                        Y=fft(S); %1 sekund signal-sekvens
                        L=length(tid(1,:));
                        P2=abs(Y/L);
                        P1=P2(1:L/2+1); %tek halve spekteret
                        P1(2:end-1)=2*P1(2:end-1); %gongar mesteparten av spekteret med 2
                        frek=Fs*(0:(L/2))/L;
                        ESD=((abs(P1)).^2)*(L/(2*Fs));
                        
                           YN=fft(SN);
                           P2N=abs(YN/L);
                           P1N=P2N(1:L/2+1);
                           P1N(2:end-1)=2*P1N(2:end-1);
                           ESDN=((abs(P1N)).^2)*(L/(2*Fs));
                           
                                   figure
                        subplot(2,1,1)
                        plot(tidsakse,data)
                        hold on
                         plot(tidsakse(tid),S,'k')
                         plot(tidsakse(tidN),SN,'r')
                         legend(['max peak=' num2str(round(20*log10(PKS/1e-6),1)) ' dB re 1 \muPa'],['1 sek seismic SEL=' num2str(round(10*log10(SEL_t/1e-12),1)) 'dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(10*log10(SEL_tN/1e-12),1)) 'dB re 1 \muPa^2\cdots'])
                      title(['1 second ' tittel], 'Interpreter', 'none')
                       xlabel('clock time, UTC')
                       ylabel('Pa')
                      set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
                       set(gca,'fontsize', 22);
                   
                        subplot(2,1,2)
%                         plot(frek,P1,'k')
%                         hold on;
%                         plot(frek,P1N,'r')
%                         legend('seismic','ambient')
%                         ylabel('Pa')
%                         xlabel('Frequency, Hz')

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
                        plot(tidsakse,data)
                        hold on
                         plot(tidsakse(tid),S,'k')
                         plot(tidsakse(tidN),SN,'r')
                         legend(['max peak=' num2str(round(20*log10(PKS/1e-6),1)) ' dB re 1 \muPa'],['1 sek seismic SEL=' num2str(round(10*log10(SEL_t/1e-12),1)) 'dB re 1 \muPa^2\cdots'],['1 sek ambient SEL=' num2str(round(10*log10(SEL_tN/1e-12),1)) 'dB re 1 \muPa^2\cdots'])
                      title(['1 second ' tittel], 'Interpreter', 'none')
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