%% starting from trimmedData
if ~exist('contrastData.mat','file')
    load('trimmedData','Cfine','perfObs','prm')
    totConc=sum(abs(Cfine),4);
    [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
    Cfine=Cfine(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
    xL = size(Cfine,1)*prm.fov(1)/1000;
    yL = size(Cfine,2)*prm.fov(2)/1000;
    zL = size(Cfine,3)*prm.fov(3)/1000;
    [Ccrs,N_el_Ccrs]=upscale(Cfine,[2 2 1 150]);
    nt=size(Ccrs,4);
    save contrastData Ccrs N_el_Ccrs  nt xL yL zL prm numIremoved numIendRemoved numJremoved numJendRemoved numKremoved numKendRemoved comprTotConc
end
