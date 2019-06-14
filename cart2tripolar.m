function [ira_tot , varargout]=cart2tripolar(im, qParams)
% Cartesian to polar pixel transform - with compact representation
% The code starts with a Cartesian image of square pixel size dx*dy=1, and
% outputs a polar image of polar pixels that has a similar amount of
% information in polar coordinates dr*dtheta~=1.  The signal in each polar
% pixel is determined via its fractional overlap with the four surrounding
% Cartesian pixels. The result is a triangular polar representation,
% because the number of polar pixels per radius per angle increase with
% radius, i.e. the larger the radius the higher the angular resolution.
% The code was originally used in the context of analyzing images with 
% symmetry around a quadrant so that functionality was kept. 
% The code support NaN valued pixels for masking.
%
% Inputs:
%  im       - an image, assumed to be centered and corrected for deformations, tilt, etc.
%             the image may contain NaN, and the function will ignore them.
%  qParams  - (optional) a vector that states the quadrants of the image to include
%             in the analysis, e.g. [2 3] will treat only the left side
%             of the image. When im is plotted using imagesc(im) then
%             qParams=1 is top right
%             qParams=2 is top left
%             qParams=3 is bottom left
%             qParams=4 is bottom right
%
% Optional Outputs:
%
% ira      - Intensity per Radius per Angle - a 2d triangular array in
%            polar coordinates of each im quadrant (pi/2 angle spread),
%            where the quadrants are represented in the 3 dimension (as 
%            pages). Area outside the triangular array is set to NaN.
% ira_tot  -  Same as ira but of the entire image that is included (up to
%            0 to 2*pi)
% spol     - taking ira_tot and stretching it to a standard 2d square array
%             using linear interpolation.
%
%
% Example:
%------------
% load('testimg.mat');
% [ira_tot ira spol]=cart2tripolar(im);
% figure(2);
% for n=1:4
%     subplot(2,2,n); imagesc(ira(:,:,n));  title(['ira - quadrant # ' num2str(n) ]);
% end
%
% figure('Position',[50 50 900 250]);
% subplot(1,3,1); imagesc(im); axis square; title('Img');
% subplot(1,3,2); imagesc(ira_tot); axis square; title('ira__tot');
% subplot(1,3,3); imagesc([0 2*pi],1:size(spol,1),spol); axis square; title('spol');
%
% __________________________________
%   Adi Natan (natan@stanford.edu)
%   Ver 1.01 , Date: Dec 6th 2017



%% Deafults
if (nargin < 2)
    qParams=1:4; % include all image quadrants unless stated otherwise.
end

if (nargin < 1) % demo mode
    qParams=1:4; % include all image quadrants unless stated otherwise.
    try load('tim.mat'); catch; im=peaks(127); end
end

im=double(im);

%% set geometry parameters
x0=ceil(size(im,1)/2); y0=ceil(size(im,2)/2); % Assuming image is centered
L = min([x0,y0]); % Length of polar image (maximal radius)
RR=(0:L);
PPR=(floor(0.5*pi*(RR+1))-1); % the # of pixels per radius for cartesian quadrant
AngleInc = (0.5*pi./PPR'); % the angle increment per radius
AngleInc(1)=0; % avoid inf at origin
%% create image quadrants
if mod(size(im,1),2)==1 %case for odd image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0:y0+L-1);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0:x0+L-1   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0:x0+L-1   ,y0:y0+L-1);
else %case for even image size
    Q(:,:,1) =  im(x0:-1:x0-L+1,y0+1:y0+L);
    Q(:,:,2) =  im(x0:-1:x0-L+1,y0:-1:y0-L+1);
    Q(:,:,3) =  im(x0+1:x0+L   ,y0:-1:y0-L+1);
    Q(:,:,4) =  im(x0+1:x0+L   ,y0+1:y0+L);
end

%%
ira = NaN(L-2,PPR(L),numel(qParams)); % initialize the  matrix

for a4n=qParams % treat each quadrant seperatly
    a4=Q(:,:,a4n);
    ira(1,1,a4n) = a4(1,1);      % origin pixel remains the same
    
    for r=2:L-2
        npr=PPR(r); % determine # polar pix in radius
        angincr=AngleInc(r);    % the angular increment per radius
        qp=0:npr;
        xp=r*sin(angincr*qp)+1;  % polar x-coordinate
        yp=r*cos(angincr*qp)+1;  % polar y-coordinate
        
        % the signal in a polar pixel at (r,qp) is determined by its
        %  fractional overlap with the four surrounding Cartesian pixels
        %  define scale fractional weight of cartesian pixels in polar pixels
        xc=round(xp);
        yc=round(yp);
        xd=1-abs(xc-xp);
        yd=1-abs(yc-yp);
        
        
        a4r=[ (a4(xc+(yc-1)*L));
            0 (a4(xc(2:npr)+(yc(2:npr)-1+(-1).^(yp(2:npr)<yc(2:npr)))*L)) 0;
            0 (a4(xc(2:npr)+(-1).^(xp(2:npr)<xc(2:npr))+yc(2:npr)*L)) 0;
            0 (a4(xc(2:npr)+(-1).^(xp(2:npr)<xc(2:npr))+L*(yc(2:npr)+(-1).^(yp(2:npr)<yc(2:npr))))) 0];
        
        c=[xd.*yd;
            0 xd(2:npr).*(1-yd(2:npr)) 0;
            0 (1-xd(2:npr)).*yd(2:npr) 0;
            0 (1-xd(2:npr)).*(1-yd(2:npr)) 0];
        
        % only when all 4 cartesian pixels are NaN the polar pixel will be NaN
        % otherwise, only the not NaN values will be considered according
        % to their fractional contribution
        
        c=bsxfun(@rdivide,c.*(~isnan(a4r)),sum(c.*(~isnan(a4r))));
        
        % gather intensity per radius per angle (ira)
        ira(r,1:npr+1,a4n) = nansum(c.*a4r);
        ira(r,all(isnan(c.*a4r)),a4n)=NaN;
    end
end

%% outputs
%  add all quadrants to one triangular array
%  (internal function)
ira_tot=addquad(ira,qParams);

% return ira and spol if asked for
nout = max(nargout,1) - 1;
if nout
    varargout{1} = ira;
    % tri-polar to square polar transform
    spol=t2s(ira,qParams);
    varargout{2} = spol;
end

if (nargin < 1) % demo mode
 figure(1);
 for n=1:4
     subplot(2,2,n); imagesc(ira(:,:,n));  title(['ira - quadrant # ' num2str(n) ]);
 end

 figure('Position',[50 50 950 260]);
 subplot(1,3,1); imagesc(im); axis square; title('Img');
 subplot(1,3,2); imagesc(ira_tot); axis square; title('ira__tot');
 subplot(1,3,3); imagesc([0 2*pi],1:size(ira,1),t2s(ira,qParams)); axis square; title('spol'); xlabel('[rad]');
end

function ira_tot=addquad(ira,qParams)
% something is buggy here with 0 replacing NaN - check and fix!
% add all quadrants to one array keeping consistent directionality
[Ms,~] = size(ira);
PPR=(floor(0.5*pi*((0:Ms)+1))-1); % calc the  # of pixels per radius in quadrant
ira_tot=NaN(size(ira,1), numel(qParams)*size(ira,2)); % preallocate

for r=Ms:-1:2
    npr=PPR(r)+1;
    clear irav
    for a4n=qParams % treat each quadrant seperatly
        if mod(a4n,2)==1 %case for odd quadrant #
            irav(:,a4n)=ira(r,1:npr,a4n);
        else %case for even quadrant #
            irav(:,a4n)=ira(r,npr:-1:1,a4n);
        end
    end
    ira_tot(r,1:npr*numel(qParams))=reshape(irav(:,qParams),1,[]);
end
ira_tot(1,1)=mean(ira(1,1,qParams),3);


function squmat=t2s(ira,qParams)
% transfer from triangular polar matrix to square polar matrix
PPR=(floor(0.5*pi*((0:size(ira,1))+1))-1); % calc the  # of pixels per radius in quadrant
%qParams=find(squeeze(~all(all(isnan(ira)))));
totPPR = numel(qParams)*(PPR+1); % total # of polar pixels per radius considering qParams

ira_tot=addquad(ira,qParams); % get the ira_tot

% create a matrix that "squares" the triangular polar matrix
squmat=zeros(size(ira_tot)); % preallocate
% the first pixel is just the intensity at the origin
squmat(1,:)=ones(1,size(ira_tot,2)).*ira_tot(1,1);
for k=2:size(ira_tot,1)
    try % in case data is all zero, then interp1 wont work
       
        
        %first map pixels in given range to new range
        smap=interp1(1:totPPR(k) , single(~isnan(ira_tot(k,1:totPPR(k)))) , linspace(1,totPPR(k),size(ira_tot,2)) ,'nearest');
        % interp only on non NaN pixels
        nr=~isnan(ira_tot(k,1:totPPR(k))); 
        squmat(k,:)=interp1(find(nr),ira_tot(k,nr),linspace(1,totPPR(k),size(ira_tot,2)),'nearest');
        squmat(k,~smap)=NaN;
       
    catch
    end
end