function SpawnSeisHydrophoneDataTreatment_testRMSogMaxogFrek(Dmeta,Tmeta,calibrationfactor)


%plukk ut tid for eksponering:
%plukkar ut liten snutt frå kvar periode
snutt=3600; %sekund
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Seismic"
        t1(teljar,:)=datetime(Tmeta(ii).Start);
        %  t2(teljar,:)=Tmeta(ii).End;
        t2(teljar,:)=t1(teljar)+datenum(0000,0,0,0,0,snutt);
        teljar=teljar+1;
    end
end

%plukk ut tid for båt kontroll (berre båt):
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Boat control"
        t1_B(teljar,:)=datetime(Tmeta(ii).Start);
        t2_B(teljar,:)=t1_B(teljar)+datenum(0000,0,0,0,0,snutt);
        %t2_B(teljar,:)=Tmeta(ii).End;
        teljar=teljar+1;
    end
end

%plukk ut tid for stille kontroll (utan båt):
teljar=1;
for ii=1:length(Tmeta)
    if string(Tmeta(ii).Treatment)=="Silent control"
        t1_S(teljar,:)=datetime(Tmeta(ii).Start);
        %t2_S(teljar,:)=Tmeta(ii).End;
        t2_S(teljar,:)=t1_S(teljar)+datenum(0000,0,0,0,0,snutt);
        teljar=teljar+1;
    end
end
%Plukk ut tider for skot som vi vil sjå nærmare på. Tider er manuelt plukka
%ut og logg i fila samples.xlsx
[~,~,Times_raw] = xlsread('samples.xlsx');


Times=cell2struct(Times_raw(2:end,:),Times_raw(1,:),2);
for ii=1:length(Times);
    
    t_p(ii,:)=datetime(Times(ii).Sample);



%% Loop over Deployments
for i=[1 2 5 6]%length(Dmeta)
    disp(Dmeta(i).Folder)
    first=1;%sjekk-parameter for å plotte rådata berre ein gong pr deployment
    k2=0;
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
            
            telj_seis=1;
            telj_sile=1;
            telj_boat=1;
            
            %For seismikk:
            if (t1(k,:)<= ref_time && t2(k,:)>= ref_time)
                
                %     tekst='seismikk';
                
                
                dat = audioread(file,'native');
                dat=detrend(double(dat));
                dat=dat*calibrationfactor(i); %calibrate
                dat=filtfilt(C,D,dat);
                
                t = ((1:length(dat))-1)*dt;
                tidsakse=ref_time+datenum(0000,0,0,0,0,t);%legg til antall sekund frå start
                
                
                if (t_p(ii,:)>=tidsakse(1) && tidsakse(1) >=t_p(ii,:)-datenum(0000,0,0,0,0,22))
                    plotSamples(dat,dt,tidsakse, ii,Dmeta(i).Folder );
                    
                end
            end
        end
    end
end
end
end

%                 [PKS,LOCS]=max(abs(dat));
%                 rms1=rms(dat);
%
%                 t_seis(telj_seis)=tidsakse(LOCS);
%                 rms_seis(telj_seis)=rms1;
%                 max_seis(telj_seis)=PKS;
%
%               %Plot ein time med skot:
% %                     if( k2==0 || k==k2)
% %                         figure(i+1000000)
% %                 plot(tidsakse,dat,'k')
% %                 hold on;
% %                 title(Dmeta(i).Folder, 'Interpreter', 'none')
% %
% %                 k2=k;
% %
% %                     first=first+1;
% %
% %                 end
%
%                    figure(i+1000)
%                 h1=plot(t_seis(telj_seis),rms_seis(telj_seis),'r*')
%                 hold on
%                                figure(i+2000)
%                 h2=plot(t_seis(telj_seis),max_seis(telj_seis),'r*')
%                 hold on
%                    figure(i+3000)
%                 h3=plot(t_seis(telj_seis),20*log10(max_seis(telj_seis)/1e-6),'r*')
%                 hold on
%
%             elseif (t1_B(k,:)<= ref_time && t2_B(k,:)>= ref_time)
%
%                     dat = audioread(file,'native');
%                 dat=detrend(double(dat));
%                 dat=dat*calibrationfactor(i); %calibrate
%                 dat=filtfilt(C,D,dat);
%                 dt = 1/ainf.SampleRate;
%                 t = ((1:length(dat))-1)*dt;
%                 tidsakse=ref_time+datenum(0000,0,0,0,0,t);%legg til antall sekund frå start
%
%                 [PKS,LOCS]=max(abs(dat));
%                 rms1=rms(dat);
%
%                 t_boat(telj_boat)=tidsakse(LOCS);
%                 rms_boat(telj_boat)=rms1;
%                 max_boat(telj_boat)=PKS;
%
%                              figure(i+1000)
%                 g1=plot(t_boat,rms_boat,'k*')
%                 hold on
%                                 figure(i+2000)
%                 g2=plot(t_boat,max_boat,'k*')
%                 hold on
%                                figure(i+3000)
%                 g2=plot(t_boat,20*log10(max_boat/1e-6),'k*')
%                 hold on
%
%
%             elseif (t1_S(k,:)<= ref_time && t2_S(k,:)>= ref_time)
%
%                       dat = audioread(file,'native');
%                 dat=detrend(double(dat));
%                 dat=dat*calibrationfactor(i); %calibrate
%                 dat=filtfilt(C,D,dat);
%                 dt = 1/ainf.SampleRate;
%                 t = ((1:length(dat))-1)*dt;
%                 tidsakse=ref_time+datenum(0000,0,0,0,0,t);%legg til antall sekund frå start
%
%                 [PKS,LOCS]=max(abs(dat));
%                 rms1=rms(dat);
%
%                 t_sile(telj_sile)=tidsakse(LOCS);
%                 rms_sile(telj_sile)=rms1;
%                 max_sile(telj_sile)=PKS;
%
%                  figure(i+1000)
%                 m1=plot(t_sile,rms_sile,'b*')
%                 hold on;
%                    figure(i+2000)
%                 m=plot(t_sile,max_sile,'b*')
%                 hold on;
%                                 figure(i+3000)
%                 m3=plot(t_sile,20*log10(max_sile/1e-6),'b*')
%                 hold on;
%
%
%
%             end
%
%
%         end
%     end
%
%              figure(i+1000)
%                 xlabel('clock time, UTC')
%                 ylabel('Pa, rms')
% %                  legend([h1 g1 m1],{'Seismic','Silent control','Boat control'})
%                 title('rms')
%                 set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
%                     'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
%                 set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
%
%                 figure(i+2000)
% %                  legend([h1 g1 m1],{'Seismic','Silent control','Boat control'})
%                 xlabel('clock time, UTC')
%                 title('maximum sound pressure')
%                 ylabel('Pa')
%                 set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
%                     'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
%                 set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
%                 title(Dmeta(i).Folder, 'Interpreter', 'none')
%
%                 figure(i+3000)
%                title('maximum sound pressure')
% %                 legend([h1 g1 m1],{'Seismic','Silent control','Boat control'})
%                 xlabel('clock time, UTC')
%                 ylabel('dB re 1\muPa')
%                 set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
%                     'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
%                 set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
%                 title(Dmeta(i).Folder, 'Interpreter', 'none')
%                 grid on
%      figure(i+1000000)
%
%     title(['Seismic exposure ' Dmeta(i).Folder])
%
%                 xlabel('clock time, UTC')
%                 ylabel('Pa')
%                 set(findobj(gcf,'type','axes'),'FontName','Calibri','FontSize',22, ...
%                     'FontWeight','Bold', 'LineWidth', 1.5,'layer','top');
%                 set(findobj(gcf, 'Type', 'Line'),'LineWidth',1.5);
%                 title(Dmeta(i).Folder, 'Interpreter', 'none')
%      grid on
% end
% end
%
%
%
%
%
%
%
