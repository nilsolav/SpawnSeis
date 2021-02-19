function SpawnSeisHydrophoneDataTreatment(Dmeta,Tmeta,calibrationfactor)


%plukk ut tid for eksponering:
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Seismic";
        t1(teljar,:)=datetime(Tmeta(ii).Start);
        t2(teljar,:)=Tmeta(ii).End;
        teljar=teljar+1;
    end
end

%plukk ut tid for båt kontroll (berre båt):
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Boat control";
        t1_B(teljar,:)=datetime(Tmeta(ii).Start);
        t2_B(teljar,:)=Tmeta(ii).End;
        teljar=teljar+1;
    end
end

%plukk ut tid for stille kontroll (utan båt):
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Silent";
        t1_S(teljar,:)=datetime(Tmeta(ii).Start);
        t2_S(teljar,:)=Tmeta(ii).End;
        teljar=teljar+1;
    end
end



%% Loop over Deployments
for i=1:6%1:length(Dmeta)
    
    % Get the file list
    filelist = dir(fullfile('./',Dmeta(i).Folder,'*wav'));
    
    for j = 1:length(filelist)
        file = fullfile(filelist(j).folder,filelist(j).name);
        ainf = audioinfo(file);
           dt = 1/ainf.SampleRate;
        %Plukker ut tid for filstart:
        
        num=str2double(regexp(ainf.Filename,'\d*','Match'));
        dato=num(end-1);
        tid=num(end);
        
        HH=floor(tid/10000);
        mm=floor((tid-(HH*10000))/100);
        ss=floor(tid-(HH*10000)-mm*100);
        yyyy=floor(dato/10000);
        SSS=0;
        MM=floor((dato -(yyyy*10000))/100);
        dd=floor(dato-(yyyy*10000)-MM*100);
        
        ref_time=datetime([num2str(dd) '-' num2str(MM) '-' num2str(yyyy) ' ' num2str(HH) ':' num2str(mm) ':' num2str(ss) '.' num2str(SSS)],'Format','dd-MM-yyyy HH:mm:ss.SSS');%tid for start av fil (antar det inntil vidare)
        
  
        
%         if i==1
%             figure(i+1000)
%             title(Dmeta(i).Folder, 'Interpreter', 'none')
%             plot(tidsakse,dat,'k')
%             hold on;
%             xlabel('Time, UTC')
%             ylabel('Pa')
%             set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
%                 'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
%             figure(i+20)
%             set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
%         end
%         
        [C,D]=butter(3,[5 5000]/(1/(2*dt)),'bandpass');
        
        %Sjekkar om snutt er innanfor eksponeringsperioden, i så fall plott peak
        for k =1:length(t1)
            
            if t1(k,:)<= ref_time && t2(k,:)>= ref_time
                
                %Les inn data for filer av interesse og kalibrerer og lagar
                %tidsakse:
                
        dat = audioread(file,'native');
        dat=detrend(double(dat));
        dat=dat*calibrationfactor(i); %calibrate
        dat=filtfilt(C,D,dat);%filtrer
        dt = 1/ainf.SampleRate;
        t = ((1:length(dat))-1)*dt;       
        tidsakse=ref_time+datenum(0000,0,0,0,0,t);%legg til antall sekund frå start
                
                figure(i+1000000)
                plot(tidsakse,dat,'k')
                hold on;
                title(Dmeta(i).Folder, 'Interpreter', 'none')
                
                
                
                
                MPH=1;%0.6*max(abs(dat));
                
                MPD=0.9; %minimum peak distance, 9.95 seconds
                l=length(dat);
                
                %  [PKS,LOCS] = findpeaks(abs(dat),l,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
                [PKS,LOCS] = findpeaks(abs(dat),'MinPeakHeight',MPH); %tar vekk sikkerhetssone
                if ~isempty(PKS)
                    if length(PKS)>1;
                        [PKS pos]=max(PKS);
                        LOCS=LOCS(pos);
                        
                    end
                    figure(i+10)
                    plot(tidsakse(LOCS),PKS,'r*')
                    hold on
                    xlabel('clock time, UTC')
                    ylabel('Pa')
                    set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                        'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
                    set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
                    grid on
                    figure(i+20)
                    
                    plot(tidsakse(LOCS),PKS,'k*')
                    hold on
                    xlabel('clock time, UTC')
                    ylabel('Pa')
                    set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                        'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
                    set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
                    title(Dmeta(i).Folder, 'Interpreter', 'none')
                    figure(i+200)
                    plot(tidsakse(LOCS),20*log10(PKS/1e-6),'k*')
                    hold on
                    xlabel('clock time, UTC')
                    ylabel('dB re 1\muPa')
                    set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                        'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
                    set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
                    title(Dmeta(i).Folder, 'Interpreter', 'none')
                    grid on
                end
            end
        end
    end
    
    
end
end




