function [allspots,roixy,framelist,params] = get_segmentation(params,datain,thr_detection)
% This function segments single emitters in the raw frames

K = params.K;
roisize = params.Mx;
maxspots = 250E3;
maxframe = 5E3;
offset = params.offset;
gain = params.gain;
maxdistance = sqrt((params.Mx^2+params.My^2))/4;

InfoImage =  imfinfo(datain);
FOV = InfoImage(1).Width;
if isfield(params,'Nframes')
    Nframes = params.Nframes;
else
    Nframes = length(InfoImage);
end
allframes = double(LoadTiff16bit(datain, [1 Nframes]));
allframes = (allframes-offset)/gain;
allframes(allframes<=0) = 1e-3;

if strcmp(params.excitation,'zstack')
    allframes_tmp = allframes;
    allframes = allframes(:,:,round(size(allframes,3)/2));
    maxdistance = sqrt((params.Mx^2+params.My^2));
    K = 1;
end
coords = zeros(maxspots,size(allframes,3)/K);
frames = zeros(maxspots,size(allframes,3)/K);

% get coordinates
fprintf(['\nPeak detection on ' num2str(size(allframes,3)) ' frames... ']); tic;
for kk = 1:size(allframes,3)/K
    
    idx = (kk-1)*K+(1:K);
    image = sum(allframes(:,:,idx),3);
    
    if strcmp(params.doetype,'none')
        ptmp = FastPeakFind(image,thr_detection,fspecial('gaussian',7,1),15,2);
    else
        ptmp = FastPeakFind(image,thr_detection,fspecial('gaussian',7,1),15,2);
    end
    
    if ~isempty(ptmp)
        
        xxtmp = round(ptmp(1:2:end));
        yytmp = round(ptmp(2:2:end));
        remove=zeros(numel(ptmp)/2,1);
        for ii = 1:numel(ptmp)/2
            distance = sqrt((xxtmp(ii)-xxtmp).^2+(yytmp(ii)-yytmp).^2);
            remove = remove|(distance<maxdistance&distance>0);
        end
        xxtmp = xxtmp(~remove);
        yytmp = yytmp(~remove);
        
        ptmp = zeros(numel(ptmp)-2*sum(remove),1);
        ptmp(1:2:end) = xxtmp;
        ptmp(2:2:end) = yytmp;
        
        numspots = numel(ptmp);
        ftmp = ones(numspots/2,1)*kk;
        
        % zero pad to match pre-allocation
        coords(:,kk) = padarray(ptmp,maxspots-numel(ptmp),'post');
        frames(:,kk) = padarray(ftmp,maxspots-numspots/2,'post');
    end
end

% turn matrix into column vector and exclude zero pads
coords = coords(coords~=0);
frames = frames(frames~=0);

ftmp = frames(:);
xtmp = round(coords(1:2:end));
ytmp = round(coords(2:2:end));
dtmp = floor(roisize/2);

fx = size(allframes,1);
fy = size(allframes,2);

ftmp = ftmp((xtmp+dtmp)<fx);
ytmp = ytmp((xtmp+dtmp)<fx);
xtmp = xtmp((xtmp+dtmp)<fx);

ftmp = ftmp((ytmp+dtmp)<fy);
xtmp = xtmp((ytmp+dtmp)<fy);
ytmp = ytmp((ytmp+dtmp)<fy);

ftmp = ftmp((xtmp-dtmp)>0);
ytmp = ytmp((xtmp-dtmp)>0);
xtmp = xtmp((xtmp-dtmp)>0);

ftmp = ftmp((ytmp-dtmp)>0);
xtmp = xtmp((ytmp-dtmp)>0);
ytmp = ytmp((ytmp-dtmp)>0);

ftmp = [ftmp(:); maxframe];

if strcmp(params.excitation,'zstack')
    allspots = zeros(roisize,roisize,size(allframes_tmp,3),numel(ftmp)-1);
else
    allspots = zeros(roisize,roisize,K,numel(ftmp)-1);
end

count = 1;
for ff = 1:size(allframes,3)/K
    idx = (ff-1)*K+(1:K);
    
    while ftmp(count) == ff
        xrange = (xtmp(count)-dtmp):(xtmp(count)+dtmp);
        yrange = (ytmp(count)-dtmp):(ytmp(count)+dtmp);
        
        if strcmp(params.excitation,'zstack')
            allspots(:,:,:,count) = allframes_tmp(yrange,xrange,:);
        else
            allspots(:,:,:,count) = allframes(yrange,xrange,idx);
        end
        
        count = count+1;
        if count>size(ftmp,1)
            break;
        end
        
    end
end

roixy = [xtmp ytmp];
framelist = ftmp(1:end-1);

% print update
fprintf(['Found ' num2str(size(allspots,4)) ' spots in ' num2str(toc,3) 's\n'])

% update parameters
params.FOV = FOV;
params.Ncfg = size(allspots,4);
params.K = size(allspots,3);
params.allalpha = repmat(params.alpha,params.Ncfg,1);

if size(allframes,3)/K == 1
    figure
    imagesc((image))
    hold on
    plot(roixy(:,1),roixy(:,2),'g+')
    pause(0.1)
end