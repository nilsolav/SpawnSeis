function SpawnSeisHydrophoneDataTreatment_testFilter(Dmeta,Tmeta,calibrationfactor)


%plukk ut tid for eksponering:
%plukkar ut liten snutt frå kvar periode
snutt=100; %sekund
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Seismic";
        t1(teljar,:)=datetime(Tmeta(ii).Start);
        %  t2(teljar,:)=Tmeta(ii).End;
        t2(teljar,:)=t1(teljar)+datenum(0000,0,0,0,0,snutt);
        teljar=teljar+1;
    end
end

%plukk ut tid for båt kontroll (berre båt):
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Boat control";
        t1_B(teljar,:)=datetime(Tmeta(ii).Start);
        t2_B(teljar,:)=t1_B(teljar)+datenum(0000,0,0,0,0,snutt);
        %t2_B(teljar,:)=Tmeta(ii).End;
        teljar=teljar+1;
    end
end

%plukk ut tid for stille kontroll (utan båt):
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Silent control";
        t1_S(teljar,:)=datetime(Tmeta(ii).Start);
        %t2_S(teljar,:)=Tmeta(ii).End;
        t2_S(teljar,:)=t1_S(teljar)+datenum(0000,0,0,0,0,snutt);
        teljar=teljar+1;
    end
end



%% Loop over Deployments
for i=2:2%length(Dmeta)
    
    % Get the file list
    filelist = dir(fullfile('./',Dmeta(i).Folder,'*wav'));
    
    for j = 1:length(filelist)
        file = fullfile(filelist(j).folder,filelist(j).name);
        ainf = audioinfo(file);
        
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
        
        
        %Sjekkar om snutt er innanfor eksponeringsperioden, i så fall plott peak
        for k =1%:length(t1)
            
            
            
            %For seismikk:
            if (t1(k,:)<= ref_time && t2(k,:)>= ref_time)|| (t1_B(k,:)<= ref_time && t2_B(k,:)>= ref_time) ||   (t1_S(k,:)<= ref_time && t2_S(k,:)>= ref_time)
                
              
                
                
                
                
                %Les inn data for filer av interesse og kalibrerer og lagar
                %tidsakse:
                
                dat = audioread(file,'native');
                dat=detrend(double(dat));
                dat=dat*calibrationfactor(i); %calibrate
                dt = 1/ainf.SampleRate;
                t = ((1:length(dat))-1)*dt;
                tidsakse=ref_time+datenum(0000,0,0,0,0,t);%legg til antall sekund frå start
                
                
                
                [C,D]=butter(3,[10 1000]/(1/(2*dt)),'bandpass');
                dat_f_10_1000=filtfilt(C,D,dat);
                
                [C,D]=butter(3,[5 5000]/(1/(2*dt)),'bandpass');
                dat_f_5_5000=filtfilt(C,D,dat);
                
                 [C,D]=butter(3,[5 1000]/(1/(2*dt)),'bandpass');
                dat_f_5_1000=filtfilt(C,D,dat);
                
                figure
                plot(tidsakse,dat)
                hold on
                plot(tidsakse,dat_f_10_1000)
                plot(tidsakse,dat_f_5_5000)
                plot(tidsakse,dat_f_5_1000)
                legend('raw','10-1000','5-30000','5-1000')
                xlabel('time, UTC')
                ylabel('Pa')
                
                
                %                 figure(i+10)
                %                 plot(tidsakse,dat,'k')
                %                 hold on;
                %                 title(Dmeta(i).Folder, 'Interpreter', 'none')
                %
                %
                %                 MPH=8;%0.6*max(abs(dat));
                %
                %                 MPD=0.9; %minimum peak distance, 9.95 seconds
                %                 l=length(dat);
                %
                %                 %  [PKS,LOCS] = findpeaks(abs(dat),l,'MinPeakHeight',MPH,'MinPeakDistance',MPD);
                %                 [PKS,LOCS] = findpeaks(abs(dat),'MinPeakHeight',MPH); %tar vekk sikkerhetssone
                %                 if ~isempty(PKS)
                %                     if length(PKS)>1;
                %                         [PKS pos]=max(PKS);
                %                         LOCS=LOCS(pos);
                [PKS,LOCS]=max(abs(dat));
                [PKS2,LOCS2]=max(abs(dat_f_10_1000));
                rms1=rms(dat);
                rms2=rms(dat_f_10_1000);
                
                
                
                figure(i+1000)
                plot(tidsakse(1),rms1,'r*')
                hold on
                plot(tidsakse(1),rms2,'b*')
                xlabel('clock time, UTC')
                ylabel('Pa, rms')
                set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                    'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
                set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
                figure(i+2000)
                
                plot(tidsakse(LOCS),PKS,'k*')
                hold on
                plot(tidsakse(LOCS2),PKS2,'b*')
                xlabel('clock time, UTC')
                ylabel('Pa')
                set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                    'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
                set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
                title(Dmeta(i).Folder, 'Interpreter', 'none')
                
                figure(i+3000)
                plot(tidsakse(LOCS),20*log10(PKS/1e-6),'k*')
                
                hold on
                plot(tidsakse(LOCS2),20*log10(PKS2/1e-6),'b*')
                xlabel('clock time, UTC')
                ylabel('dB re 1\muPa')
                set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
                    'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
                set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
                title(Dmeta(i).Folder, 'Interpreter', 'none')
            end
        end
    end
end
end







