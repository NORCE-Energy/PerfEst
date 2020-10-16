function [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,activeBlocksCompressed]=removeUnactiveLayers(activeBlocks)

% remove unactive layers
numIremoved=0;
while sum(sum(activeBlocks(1,:,:),'omitnan'))==0 % && sum(sum(bottomBlocks(1,:,:)))==0
    activeBlocks(1,:,:)=[];
    %bottomBlocks(1,:,:)=[];
    numIremoved=numIremoved+1;
end

numIendRemoved=0;
while sum(sum(activeBlocks(end,:,:),'omitnan'))==0 % && sum(sum(bottomBlocks(end,:,:)))==0
    activeBlocks(end,:,:)=[];
    %bottomBlocks(end,:,:)=[];
    numIendRemoved=numIendRemoved+1;
end

numJremoved=0;
while sum(sum(activeBlocks(:,1,:),'omitnan'))==0 % && sum(sum(bottomBlocks(:,1,:)))==0
    activeBlocks(:,1,:)=[];
    %bottomBlocks(:,1,:)=[];
    numJremoved=numJremoved+1;
end

numJendRemoved=0;
while sum(sum(activeBlocks(:,end,:),'omitnan'))==0 % && sum(sum(bottomBlocks(:,end,:)))==0
    activeBlocks(:,end,:)=[];
    %bottomBlocks(:,end,:)=[];
    numJendRemoved=numJendRemoved+1;
end

numKremoved=0;
while sum(sum(activeBlocks(:,:,1),'omitnan'))==0
    activeBlocks(:,:,1)=[];
    numKremoved=numKremoved+1;
end

numKendRemoved=0;
while sum(sum(activeBlocks(:,:,end),'omitnan'))==0
    activeBlocks(:,:,end)=[];
    numKendRemoved=numKendRemoved+1;
end
activeBlocksCompressed=activeBlocks;