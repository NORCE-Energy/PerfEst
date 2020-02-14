%% original size
size(contrast)

%% upscaled contrast
upscaledContrast=zeros(110,120,30,90);
maxX=size(contrast,1);
maxY=size(contrast,2);
maxZ=size(contrast,3);

nx=ceil(size(contrast,1)/size(upscaledContrast,1))
ny=ceil(size(contrast,2)/size(upscaledContrast,2))
nz=ceil(size(contrast,3)/size(upscaledContrast,3))

for I=1:size(upscaledContrast,4)
    for J=1:size(upscaledContrast,3) % start from 2
        for K=1:size(upscaledContrast,2) 
            for L=1:size(upscaledContrast,1)
                upperX=min([L*nx,maxX]);upperY=min([K*ny,maxY]);upperZ=min([1+J*nz,maxZ]);
                upscaledContrast(L,K,J,I)=sum(sum(sum(contrast(1+nx*(L-1):upperX,1+ny*(K-1):upperY,2+nz*(J-1):upperZ ,I))));
            end
        end
    end
    I
end
            
activeBlocks=zeros(size(upscaledContrast(:,:,:,1)));

for I=1:size(upscaledContrast,1)
    for J=1:size(upscaledContrast,2)
        for K=1:size(upscaledContrast,3)
            activeBlocks(I,J,K)=max(upscaledContrast(I,J,K,:));
        end
    end
end
activeBlocks(activeBlocks>0)=1;

%% upscale bottom layer
bottomContrast=zeros(size(upscaledContrast,1),size(upscaledContrast,2),1,90);
maxX=size(contrast,1);
maxY=size(contrast,2);
maxZ=1;

nx=ceil(size(contrast,1)/size(bottomContrast,1))
ny=ceil(size(contrast,2)/size(bottomContrast,2))
multFactor=nz
%nz=1; use the previous one as a multiplication factor
for I=1:size(bottomContrast,4)
    for J=1:size(bottomContrast,3) %
        for K=1:size(bottomContrast,2) 
            for L=1:size(bottomContrast,1)
                upperX=min([L*nx,maxX]);upperY=min([K*ny,maxY]);upperZ=min([1+J,maxZ]);
                bottomContrast(L,K,J,I)=multFactor*sum(sum(sum(contrast(1+nx*(L-1):upperX,1+ny*(K-1):upperY,1+nz*(J-1):upperZ ,I))));
            end
        end
    end
    I
end
            
bottomBlocks=zeros(size(bottomContrast(:,:,:,1)));
for I=1:size(bottomContrast,1)
    for J=1:size(bottomContrast,2)
        for K=1:size(bottomContrast,3)
            bottomBlocks(I,J,K)=max(bottomContrast(I,J,K,:));
        end
    end
end
bottomBlocks(bottomBlocks>0)=1;

%% remove unactive layers
numIremoved=0;
while sum(sum(activeBlocks(1,:,:)))==0 && sum(sum(bottomBlocks(1,:,:)))==0
    activeBlocks(1,:,:)=[];
    bottomBlocks(1,:,:)=[];
    numIremoved=numIremoved+1;
end

numIendRemoved=0;
while sum(sum(activeBlocks(end,:,:)))==0 && sum(sum(bottomBlocks(end,:,:)))==0
    activeBlocks(end,:,:)=[];
    bottomBlocks(end,:,:)=[];
    numIendRemoved=numIendRemoved+1;
end

numJremoved=0;
while sum(sum(activeBlocks(:,1,:)))==0 && sum(sum(bottomBlocks(:,1,:)))==0
    activeBlocks(:,1,:)=[];
    bottomBlocks(:,1,:)=[];
    numJremoved=numJremoved+1;
end

numJendRemoved=0;
while sum(sum(activeBlocks(:,end,:)))==0 && sum(sum(bottomBlocks(:,end,:)))==0
    activeBlocks(:,end,:)=[];
    bottomBlocks(:,end,:)=[];
    numJendRemoved=numJendRemoved+1;
end

numKremoved=0;
while sum(sum(activeBlocks(:,:,1)))==0
    activeBlocks(:,:,1)=[];
    numKremoved=numKremoved+1;
end

numKendRemoved=0;
while sum(sum(activeBlocks(:,:,end)))==0
    activeBlocks(:,:,end)=[];
    numKendRemoved=numKendRemoved+1;
end

%% compress upscaled contrast
upscaledContrast=upscaledContrast(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
upscaledContrast(:,:,2:end+1,:)=upscaledContrast;
upscaledContrast(:,:,1,:)=bottomContrast(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,1,:);
activeBlocks(:,:,2:end+1)=activeBlocks;
activeBlocks(:,:,1)=bottomBlocks;


