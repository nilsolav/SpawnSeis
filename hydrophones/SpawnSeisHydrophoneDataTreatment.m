function SpawnSeisHydrophoneDataTreatment(Dmeta,Tmeta)
%% Loop over Deployments
for i=2%1:length(Dmeta)
    
    % Get the file list
    filelist = dir(fullfile('./',Dmeta(i).Folder,'*wav'));
    disp(filelist)

    for j = 2%length(filelist)
        file = fullfile(filelist(j).folder,filelist(j).name);
        ainf = audioinfo(file);
        dat = audioread(file);
        dt = 1/ainf.SampleRate;
        t = ((1:length(dat))-1)*dt;
        plot(t,dat)
    end
end



