function [DEM,ix] = largestinscribedgrid(DEM,varargin)

%LARGESTINSCRIBEDGRID Find and crop the largest grid with no nans
%
% Syntax
%
%     DEMc = largestinscribedgrid(DEM)
%     DEMc = largestinscribedgrid(DEM,pn,pv,...)
%     [DEMc,ix] = ...
%
% Description
%
%     This function finds the grid with the maximum area which contains no
%     nans. For example, projection of a grid into another coordinate
%     reference system will result in nans along the grid borders. To
%     remove the nans, you will need to identify the grid that is as large
%     as possible but does not contain nans. largestinscribedgrid uses the
%     function LargestRectangle by Peter Seibold (2020).
%
% Input arguments
%
%     DEM     GRIDobj with nans
%
%     Parameter name/value pairs
%
%     'holes'  {false} or true. If true, then the resulting may contain
%              holes of nans (by chance, these may actually be on the 
%              border of the cropped grid).
%     
%   
% Output arguments
%
%     DEMc    Cropped DEM with no nans
%     ix      linear indices of corners (can be used to crop other grids)
%
% Example
%
%     DEM = readexample('taiwan');
%     [DEMc,ix] = largestinscribedgrid(DEM);
%     imagesc(DEM)
%     hold on
%     imageschs(DEMc)
%     [x,y] = ind2coord(DEM,ix);
%     plot(x,y,'+r')
%
% Reference: Peter Seibold (2020). Largest inscribed rectangle square or 
%            circle, MATLAB Central File Exchange. 
%
% See also: GRIDobj/crop, readopentopo, reproject2utm
%
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 14. August, 2023

p = inputParser;
addParameter(p,'holes',false)
parse(p,varargin{:})

% Find non-nan region
I = ~isnan(DEM.Z);
% Remove holes if required
if p.Results.holes
    I = imfill(I,'holes');
end

% Find largest rectangle without rotation
L = LargestRectangle(I,90,0,0,89.999,0);
% Convert subscripts to indices
ix = sub2ind(DEM.size,L(2:end,2),L(2:end,1));
% And clip the DEM
DEM = crop(DEM,ix);


end

function LRout=LargestRectangle(image,varargin)
%Function to find the largest inscribed rectangle in an arbitrary shape with multiple holes.
%INPUT:
% image:        Image, RGB, grey or BW. By preference BW.
% RotationStep: Default: 5°
%               Rotation Step in degrees. In order to find tilted rectangles, the image is rotated.
%               Range: 0 < RotationStep <= 90. 
%               If RotationStep>(LastAngle-FirstAngle) then the image is rotated once.
% iterate:      Default: 1 (with iteration)
%               If 0 then no iteration, if 1 then iteration. 
%               No iteration if (FirstAngle == LastAngle)
%               It might not find the largest rectangle in the range determined by the RotationStep,
%               it might lock on any local maximum within a rotation step.
%FirstAngle:    Default: 0°
%               First angle:  0 <= FirstAngle < 90
%LastAngle:     Default: 89.9999°
%               Last angle: FirstAngle <= LastAngle <= 90 
%Graphic:       Default: 1 (Plot graphic), 0: no graphic
%OUTPUT LRout:
% 1st row: Area of the largest rectangle in px, Rotation angle in degrees counterclockwise
% 2nd row: x,y of top corner of largest rectangle
% 3rd row: x,y of right corner of largest rectangle
% 4th row: x,y of bottom corner of largest rectangle
% 5th row:x,y of left corner of largest rectangle
%EXAMPLES:
% LRout=LargestRectangle(myImage)% Run with default values
% LRout=LargestRectangle(myImage,0.1,0)% Run with small rotation steps, no iteration
% LRout=LargestRectangle(myImage,0,0,0,0)% For axis parallel rectangle
% Example for rotation step=5, no iteration, starting at 10°, ending at 20°, no graphic:
%    LRout=LargestRectangle(myImage,5,0,10,20,0)
%REMARK:
% Any hole (black pixel), even only one pixel wide, 
%  will not be inside the largest rectangle
%Revisions:
%19.05.20: removed bug when left or top border of input image has a white pixel
inArg=[5,1,0,89.9999,1];%Default values, you may want to change it
inArg(1:nargin-1)=[varargin{1:nargin-1}];
RotationStep=inArg(1);
iterate=inArg(2);
StartAngle=inArg(3);
EndAngle=inArg(4);
numRotSteps=floor((EndAngle-StartAngle)/RotationStep)+1;
if EndAngle==StartAngle
  iterate=0;% NO iteration
  RotationStep=1;% dummy value >0
  numRotSteps=1;
end
if iterate && floor((EndAngle-StartAngle)/RotationStep)+1<10
  % make rotatationstep smaller and keep user selected values
    RotationStep=RotationStep/ceil(10*RotationStep/(EndAngle-StartAngle));
    numRotSteps=floor((EndAngle-StartAngle)/RotationStep)+1;
end
if ~islogical(image)
  image=im2bw(image);
end
ImBWc=image;
%add 1px black border arround ImBWc:
%You gain approx. 1% area if you remove 'logical', but 3 times slower
ImBWc=logical([zeros(1,size(ImBWc,2)+2);...%add top border
  zeros(size(ImBWc,1),1),ImBWc,zeros(size(ImBWc,1),1);...%add left+right border  + image
  zeros(1,size(ImBWc,2)+2)]);%add bottom border
[rows, columns] = find(ImBWc);
if isempty(rows)
  %totaly black image
  disp ('Totaly black image! No rectangle found.');
  LRout=[0,0;1,1;1,1;1,1;1,1];
  return;
end
%determine smallest black border arround image
leftCrop = min(columns)-1;
topCrop = min(rows)-1;
rightCrop = max(columns)+1;
bottomCrop = max(rows)+1;
%crop image:
ImBWc=ImBWc(topCrop:bottomCrop,leftCrop:rightCrop);
%correct origin for later use (1px was added as border):
leftCrop = leftCrop-1;
topCrop = topCrop-1;
LR=[0,0,0,0,0,0,0,0,0];%area,xmin,xmax,ymin,ymax, angle, sizeImR(1), sizeImR(2),AreaNow
%Main loop, rotate image get border, get largest rectangle
AN=zeros(numRotSteps,2);% Angle, LR, used for iterate only
ANi=0;
for angle=StartAngle:RotationStep:EndAngle
    [border,ColLims,RowLims,sizeImR]=GetContour(ImBWc,angle);
    if border(1)==0
      %Only 1 px areas
      [border,ColLims,RowLims,sizeImR]=GetContour(ImBWc,0);
      LR=getLR(LR,border,ColLims,RowLims);
      LR(6:8)=[0,sizeImR];
      iterate=0;
      break;%Exit for loop
    end    
    areaOld=LR(1);
    LR=getLR(LR,border,ColLims,RowLims);
    ANi=ANi+1;
    AN(ANi,:)=[angle,LR(9)];    
    if LR(1)>areaOld
        LR(6:8)=[angle,sizeImR];
    end
end
if iterate 
  %iterate
  MinRotationStep=min(0.01,RotationStep/10);
  AN=[AN;AN(end,1)+RotationStep,0];
  ANi=ANi+1;
  LimitDiff=0.05*0.85^-log2(RotationStep/15);
  angles=[AN(1,1),AN(end,1)];
  Transit=true;%for first iteration at <0.5°
  while RotationStep>MinRotationStep 
    RotationStep2=RotationStep/2;
    numelAN=size(AN,1);
    LimitDiff=LimitDiff*0.85;%decreasing last number results in narrower search range 
    Limit=LR(1)*(1-LimitDiff);    
    iLold=find(AN(:,1)>angles(1)-RotationStep2,1,'first');
    iRold=find(AN(:,1)>angles(end)-RotationStep2,1,'first');
    iLRmaxM=find(AN(iLold:iRold,2)==LR(1))+iLold-1;%multiple iLRmax
    %check if the LRmax are at around 0° AND around 90° (border problem!)
    if AN(iLRmaxM(1,1))<15 && AN(iLRmaxM(end,1))>85
      %border problem!
      iLRmaxML=iLRmaxM(AN(iLRmaxM(:),1)<45);
      iLRmaxMR=iLRmaxM(iLRmaxML+1:end);
      %take the side with the most max LR
      if numel(iLRmaxML)>=numel(iLRmaxMR)
        iLRmaxM=iLRmaxML;
      else
        iLRmaxM=iLRmaxMR;
      end
    end
    iLRmax=iLRmaxM(floor((numel(iLRmaxM)+1)/2));%take position of midle LRmax
    if Transit && RotationStep<0.5
      %correct iLold and iRold
      if inArg(1)<0.5
        %no previous iteration
        iLstart=max(1,iLRmaxM(1)-20);
        iRend=min(numelAN,iLRmaxM(end)+20);
        iLold=find(AN(iLstart:iLRmaxM(1),2)<Limit,1,'last')+iLstart-1;%Left index
        iRold=find(AN(iLRmaxM(end):iRend,2)<Limit,1,'first')+iLRmaxM(end)-1;%Right index
        if isempty(iLold);iLold=iLstart;end
        if iLold<1; iLold=1;end
        if isempty(iRold);iRold=iRend;end
        if iRold>iRend;iRold=iRend;end        
      else
        iLold=find(AN(:,1)==angles(find(angles>AN(iLRmaxM(1),1)-5,1,'first')));
        iRold=find(AN(:,1)==angles(find(angles<AN(iLRmaxM(end),1)+5,1,'last')));
      end
      Transit=false;
    end
    iLstart=max(iLRmax-25,iLold);
    iRend=min(iLRmax+25,iRold); 
    if RotationStep<0.5 %last number by experiment
      %include limited range arround max peak and second largest peak
      ANs=sort(AN(iLstart:iRend,2),'descend');
      ANs=ANs(ANs<LR(1));%exclude (multiple) LR(1)
      if isempty(ANs)
       LimitE=LR(1)-1;
      else
        LimitE=ANs(1)-1;%second largest-1
      end
      iLE=find(AN(iLstart:iLRmax-1,2)>LimitE,1,'first')+iLstart-2;
      if iLE<1;iLE=iLstart;end
      iRE=find(AN(iLRmax+1:iRend,2)>LimitE,1,'last')+iLRmax+1;%Right index 
      if iRE>numelAN;iRE=numelAN;end
      iL=find(AN(iLstart:iLRmax,2)<Limit,1,'last')+iLstart-1;%Left index
      if isempty(iL) || iL<1;iL=iLstart;end
      iR=find(AN(iLRmax:iRend,2)<Limit,1,'first')+iLRmax-1;%Right index 
      if isempty(iR) ||iR>iRend;iR=iRend;end
      if iLE<iL
        iL=iLE;
      end
      if iRE>iR
        iR=iRE;
      end
      maxRotIt=min(22,round(6.5+1.7./sqrt(RotationStep)));%9 at 0.4°, 22 at 0.1°
      maxRot=maxRotIt;
        while iR-iL>maxRot
            %delete smaller LR
            if AN(iL,2)<= AN(iR,2) && iL<iLRmaxM(1)
              iL=iL+1;
            elseif iR>iLRmaxM(end)
              iR=iR-1;
            else
              iL=max(iL,iLRmax-12);
              iR=min(iR,iLRmax+12);
              break;
            end
            if iL==iLRmaxM(1) ||  iR==iLRmaxM(end)
               %propably unsymetric peak position
              maxRot=24;
            end
        end      
      if iR-iL<4 || (iL>=iLRmaxM(1) && iR<=iLRmaxM(end))
        %only a few rot steps, increase by 2
        angles=AN(iL,1)-RotationStep2:RotationStep:AN(iR,1)+1.1*RotationStep2;      
      elseif iL>=iLRmaxM(1)
        angles=AN(iL,1)-RotationStep2:RotationStep:AN(iR,1);
      elseif iR<=iLRmaxM(end)
        angles=AN(iL,1)+RotationStep2:RotationStep:AN(iR,1)+1.1*RotationStep2;
      else
        %regular case
        angles=AN(iL,1)+RotationStep2:RotationStep:AN(iR,1);
      end      
    else
      %RotationStep>=0.5, select largest peaks
      iL=find(AN(iLstart:iLRmax,2)>Limit,1,'first')+iLstart-2;%Left index
      iR=find(AN(iLRmax:iRend,2)>Limit,1,'last')+iLRmax;%Right index
      if isempty(iL);iL=iLRmax-2;end
      if iL<1; iL=1;end
      if isempty(iR);iR=iLRmax+2;end
      if iR>iRend;iR=iRend;end
      LimitLow=2*Limit-LR(1);
      angles=AN(iL:iR,:);
      angles=angles(angles(:,2)>LimitLow,:);%exclude small LRs      
      maxRot=7;
      if length(angles)>maxRot
        angles=sortrows(angles,-2);
        angles=angles(1:maxRot,1);
      else
        angles=angles(:,1);
      end      
      angles=sort([angles-RotationStep2;angles+RotationStep2]);
      %exclude duplicate angles:
      angles=[angles(abs(angles(1:end-1)-angles(2:end))>RotationStep2);angles(end)];      
    end% if RotationStep<0.5
    for j=1:numel(angles)
        angle=angles(j);
        [border,ColLims,RowLims,sizeImR]=GetContour(ImBWc,angle);
        areaOld=LR(1);
        LR=getLR(LR,border,ColLims,RowLims);
        ANi=ANi+1;
        AN(ANi,:)=[angle,LR(9)];    
        if LR(1)>areaOld
            LR(6:8)=[angle,sizeImR];
        end
    end
    AN=sortrows(AN); 
    RotationStep=RotationStep2;
  end
end%End iterate
%Prepaire LRout
ImRcenter=([LR(7),LR(8)]+1)/2;
ImBWcCenter=(size(ImBWc)+1)/2;
xy=[[LR(2),LR(3),LR(3),LR(2)]-ImRcenter(2);[LR(4),LR(4),LR(5),LR(5)]-ImRcenter(1)];
phi=LR(6)/180*pi;
cosPhi=cos(phi);
sinPhi=sin(phi);
RM= [cosPhi -sinPhi; sinPhi cosPhi];
xy=RM*xy;
x=xy(1,:)+ImBWcCenter(2)+leftCrop-1;
y=xy(2,:)+ImBWcCenter(1)+topCrop-1;
LRout=[LR(1),LR(6);x(1),y(1);x(2),y(2);x(3),y(3);x(4),y(4)];
if inArg(5)
  plotResult(image,LR,x,y);
end

function [borderD2,ColLims,RowLims,sizeImR]=GetContour(ImBWc,angle)
  ImBWc=imrotate(ImBWc,angle,'bilinear');
  sizeImR=size(ImBWc);
  sizeImR1=sizeImR(1);
  boundary = bwperim(ImBWc);
  [row,col] = find(boundary);
  borderLinIdxL=(col-2)*sizeImR1+row;% indices for left side of border point 
  borderLinIdxR=borderLinIdxL+2*sizeImR1;% indices for right of border points
  borderLinIdxT=borderLinIdxL+sizeImR1-1;% indices for above of border points 
  borderLinIdxB=borderLinIdxT+2;% indices for below of border points
  %top-bottom limits:
  borderLinIdx=[borderLinIdxT;borderLinIdxB];
  indx=unique(borderLinIdx(ImBWc(borderLinIdx)==0));
  c=floor((indx-1)/sizeImR1)+1;
  borderD=[c,indx-(c-1)*sizeImR1];%x,y
  [~,mm]=mode(borderD(:,1));
   numRows=borderD(end,1);
   indx=borderD(:,1); 
   indx3=indx;
   ColLims=NaN(numRows,mm);
   for i=2:numel(indx)
     if indx(i)==indx(i-1)
        indx3(i)=indx3(i-1)+numRows;
     else
       indx3(i)=indx(i);
     end
   end
   ColLims(indx3)=borderD(:,2);   
  %Left-right limits
  borderLinIdx=[borderLinIdxL;borderLinIdxR];
  indx=unique(borderLinIdx(ImBWc(borderLinIdx)==0));
  c=floor((indx-1)/sizeImR1)+1;
  borderD=sortrows([c,indx-(c-1)*sizeImR1],2);%x,y
  [~,mm]=mode(borderD(:,2));
   numRows=borderD(end,2);
   indx=borderD(:,2); 
   indx3=indx;
   RowLims=NaN(numRows,mm);
   for i=2:numel(indx)
     if indx(i)==indx(i-1)
        indx3(i)=indx3(i-1)+numRows;
     else
       indx3(i)=indx(i);
     end
   end
   RowLims(indx3)=borderD(:,1);
  %find top border px only
  indxT=borderLinIdxT(ImBWc(borderLinIdxT)==0);
  %find left border px only
  indxL=borderLinIdxL(ImBWc(borderLinIdxL)==0);
  if numel(indxT)<numel(indxL)
    %Use top border only for evaluation
    c=floor((indxT-1)/sizeImR1)+1;
    borderD2=[c,indxT-(c-1)*sizeImR1+1];
  else
    c=floor((indxL-1)/sizeImR1)+1;
    borderD2=[c+1,indxL-(c-1)*sizeImR1];
  end
end

function LR=getLR(LR,border,ColLims,RowLims)
  %get largest rectangle
  AreaNow=0;%only for iteration
  for i=1:size(border,1)
    %loop border points
    xNow=border(i,1);
    yNow=border(i,2);
    Indx=1;while ColLims(xNow,Indx)<yNow;Indx=Indx+1;end%much faster than find.m    
    yLimit1=ColLims(xNow,Indx-1)+1;yLimit2=ColLims(xNow,Indx)-1;%ylimits @xNow
    Indx=1;while RowLims(yNow,Indx)<xNow;Indx=Indx+1;end    
    xLimit1=RowLims(yNow,Indx-1)+1;xLimit2=RowLims(yNow,Indx)-1;%xlimits @yNow
    heightMax=yLimit2-yLimit1+1;
    widthMax=xLimit2-xLimit1+1;
    if widthMax*heightMax>LR(1)
      if (yNow-yLimit1+1)*(yLimit2-yNow+1)<(xNow-xLimit1+1)*(xLimit2-xNow+1) 
        %Max Number of Operation y direction smaller than for x direction 
        if yNow==yLimit1
          %Border pixel at top  ===========================================
            xLimitTempT1=xLimit1;
            xLimitTempT2=xLimit2;
            xLimitTempB1=xLimit1;
            xLimitTempB2=xLimit2;
            Indx=1;while RowLims(yNow,Indx)<xNow;Indx=Indx+1;end        
            xT2=RowLims(yNow,Indx)-1;xT1=RowLims(yNow,Indx-1)+1;%faster in one line
            if xT1>xLimitTempT1;xLimitTempT1=xT1;end
            if xT2<xLimitTempT2;xLimitTempT2=xT2;end
            if (xLimitTempT2-xLimitTempT1+1)*heightMax<LR(1)  
              break;
            end    
            for yB=yNow:yLimit2
              Indx=1;while RowLims(yB,Indx)<xNow;Indx=Indx+1;end 
              xB2=RowLims(yB,Indx)-1;xB1=RowLims(yB,Indx-1)+1;%faster in one line
              if xB1>xLimitTempB1;xLimitTempB1=xB1;end
              if xB2<xLimitTempB2;xLimitTempB2=xB2;end          
              if xLimitTempB1>xLimitTempT1;xL=xLimitTempB1;else xL=xLimitTempT1;end
              if xLimitTempB2>xLimitTempT2;xR=xLimitTempT2;else xR=xLimitTempB2;end              
              if (xR-xL+1)*heightMax<LR(1)
                break;
              end
              area=(xR-xL+1)*(yB-yNow+1);
              if area>AreaNow
                AreaNow=area;
                if area>LR(1)
                  LR(1:5)=[area,xL,xR,yNow,yB];
                end                
              end               
            end%for yB=yNow:yLimit2        
      elseif yNow==yLimit2
          %Border pixel at bottom  ========================================
          xLimitTempB1=xLimit1;
          xLimitTempB2=xLimit2;
          xLimitTempT1=xLimit1;
          xLimitTempT2=xLimit2;
          Indx=1;while RowLims(yNow,Indx)<xNow;Indx=Indx+1;end    
          xB2=RowLims(yNow,Indx)-1;
          xB1=RowLims(yNow,Indx-1)+1;
          if xB1>xLimitTempB1;xLimitTempB1=xB1;end
          if xB2<xLimitTempB2;xLimitTempB2=xB2;end
          if (xLimitTempB2-xLimitTempB1+1)*heightMax>LR(1)
            for yT=yNow:-1:yLimit1
              Indx=1;while RowLims(yT,Indx)<xNow;Indx=Indx+1;end
              xT2=RowLims(yT,Indx)-1;
              xT1=RowLims(yT,Indx-1)+1;        
              if xT1>xLimitTempT1;xLimitTempT1=xT1;end
              if xT2<xLimitTempT2;xLimitTempT2=xT2;end 
              if xLimitTempB1>xLimitTempT1;xL=xLimitTempB1;else xL=xLimitTempT1;end
              if xLimitTempB2<xLimitTempT2;xR=xLimitTempB2;else xR=xLimitTempT2;end
              if (xR-xL+1)*heightMax<LR(1)
                break;
              end
              area=(xR-xL+1)*(yNow-yT+1);
              if area>AreaNow
                AreaNow=area;
                if area>LR(1)
                  LR(1:5)=[area,xL,xR,yT,yNow];
                end
              end            
            end%
          end%if (xLimitTempB2-xLimitTempB1+1)*heightMax>LR(1)
        else
          %Border pixel between top and bottom ============================
          xLimitTempT1=xLimit1;
          xLimitTempT2=xLimit2;
          %prepare data for inner loop
          xB1a=zeros(yLimit2-yNow+1,1);
          xB2a=xB1a;
          for yB=yNow:yLimit2
            Indx=1;while RowLims(yB,Indx)<xNow;Indx=Indx+1;end 
            xB2a(yB-yNow+1)=RowLims(yB,Indx)-1;xB1a(yB-yNow+1)=RowLims(yB,Indx-1)+1;%faster in one line            
          end  
          for yT=yNow:-1:yLimit1
            xLimitTempB1=xLimit1;
            xLimitTempB2=xLimit2;
            Indx=1;while RowLims(yT,Indx)<xNow;Indx=Indx+1;end        
            xT2=RowLims(yT,Indx)-1;xT1=RowLims(yT,Indx-1)+1;%faster in one line
            if xT1>xLimitTempT1;xLimitTempT1=xT1;end
            if xT2<xLimitTempT2;xLimitTempT2=xT2;end
            if (xLimitTempT2-xLimitTempT1+1)*heightMax<LR(1)  
              break;
            end 
            for yB=yNow:yLimit2
              if xB1a(yB-yNow+1)>xLimitTempB1;xLimitTempB1=xB1a(yB-yNow+1);end
              if xB2a(yB-yNow+1)<xLimitTempB2;xLimitTempB2=xB2a(yB-yNow+1);end          
              if xLimitTempB1>xLimitTempT1;xL=xLimitTempB1;else xL=xLimitTempT1;end
              if xLimitTempB2>xLimitTempT2;xR=xLimitTempT2;else xR=xLimitTempB2;end              
              if (xR-xL+1)*(yLimit2-yT+1)<LR(1)
                break;
              end              
              area=(xR-xL+1)*(yB-yT+1);
              if area>AreaNow
                AreaNow=area;
                if area>LR(1)
                  LR(1:5)=[area,xL,xR,yT,yB];
                end
              end               
            end%for yB=yNow:yLimit2
          end%for yT=yNow:-1:yLimit1           
        end%
      else%
        %Max Number of Operation x direction smaller than for y direction ==================
        if xNow==xLimit1
          %Border pixel on left side ======================================
            yLimitTempL1=yLimit1;%ylimit Left
            yLimitTempL2=yLimit2;%ylimit Left
            yLimitTempR1=yLimit1;
            yLimitTempR2=yLimit2;
            Indx=1;while ColLims(xNow,Indx)<yNow;Indx=Indx+1;end     
            yL2=ColLims(xNow,Indx)-1;yL1=ColLims(xNow,Indx-1)+1;%faster in one line
            if yL1>yLimitTempL1;yLimitTempL1=yL1;end
            if yL2<yLimitTempL2;yLimitTempL2=yL2;end
            if widthMax*(yLimitTempL2-yLimitTempL1+1)>LR(1)
              for xR=xNow:xLimit2
                Indx=1;while ColLims(xR,Indx)<yNow;Indx=Indx+1;end              
                yR2=ColLims(xR,Indx)-1;yR1=ColLims(xR,Indx-1)+1;%faster in one line
                if yR1>yLimitTempR1;yLimitTempR1=yR1;end%much faster than yR1=max([yLimitTempR1,yR1]);
                if yR2<yLimitTempR2;yLimitTempR2=yR2;end            
                if yLimitTempR1>yLimitTempL1;yT=yLimitTempR1;else yT=yLimitTempL1;end
                if yLimitTempR2>yLimitTempL2;yB=yLimitTempL2;else yB=yLimitTempR2;end  
                if widthMax*(yB-yT+1)<LR(1)
                  break;
                end
                area=(xR-xNow+1)*(yB-yT+1);
                if area>AreaNow
                  AreaNow=area;
                  if area>LR(1)
                    LR(1:5)=[area,xNow,xR,yT,yB];
                  end                  
                end                     
              end%
            end%if widthMax*(yLimitTempL2-yLimitTempL1+1)>LR(1)
         elseif xNow==xLimit2%if xNow==xLimit1
          %Border pixel on right side =====================================
          yLimitTempR1=yLimit1; 
          yLimitTempR2=yLimit2; 
            Indx=1;while ColLims(xNow,Indx)<yNow;Indx=Indx+1;end              
            yR2=ColLims(xNow,Indx)-1;yR1=ColLims(xNow,Indx-1)+1;%faster in one line
            if yR1>yLimitTempR1;yLimitTempR1=yR1;end
            if yR2<yLimitTempR2;yLimitTempR2=yR2;end             
            if widthMax*(yLimitTempR2-yLimitTempR1+1)<LR(1)  
              break;
            end
            yLimitTempL1=yLimit1;
            yLimitTempL2=yLimit2;
            for xL=xNow:-1:xLimit1  
              Indx=1;while ColLims(xL,Indx)<yNow;Indx=Indx+1;end
              yL2=ColLims(xL,Indx)-1;yL1=ColLims(xL,Indx-1)+1;%faster in one line
              if yL1>yLimitTempL1;yLimitTempL1=yL1;end
              if yL2<yLimitTempL2;yLimitTempL2=yL2;end
              if yLimitTempR1>yLimitTempL1;yT=yLimitTempR1;else yT=yLimitTempL1;end
              if yLimitTempR2>yLimitTempL2; yB=yLimitTempL2;else yB=yLimitTempR2;end 
              if widthMax*(yB-yT+1)<LR(1)
                break;
              end
              area=(xNow-xL+1)*(yB-yT+1);
              if area>AreaNow
                AreaNow=area;
                if area>LR(1)
                  LR(1:5)=[area,xL,xNow,yT,yB];
                end
              end             
            end%for xL=xNow:-1:xLimit1
        else
          %Border pixel between left and right side =======================
          yLimitTempR1=yLimit1; 
          yLimitTempR2=yLimit2;
          %prepare data for inner loop
          yL1a=zeros(xNow-xLimit1+1,1);
          yL2a=yL1a;
          for xL=xNow:-1:xLimit1
              Indx=1;while ColLims(xL,Indx)<yNow;Indx=Indx+1;end
              yL2a(xL-xLimit1+1)=ColLims(xL,Indx)-1;yL1a(xL-xLimit1+1)=ColLims(xL,Indx-1)+1;            
          end
          for xR=xNow:xLimit2
            Indx=1;while ColLims(xR,Indx)<yNow;Indx=Indx+1;end              
            yR2=ColLims(xR,Indx)-1;yR1=ColLims(xR,Indx-1)+1;%faster in one line
            
            if yR1>yLimitTempR1;yLimitTempR1=yR1;end
            if yR2<yLimitTempR2;yLimitTempR2=yR2;end             
            if widthMax*(yLimitTempR2-yLimitTempR1+1)<LR(1)  
              break;
            end
            yLimitTempL1=yLimit1;
            yLimitTempL2=yLimit2;
            for xL=xNow:-1:xLimit1  
              if yL1a(xL-xLimit1+1)>yLimitTempL1;yLimitTempL1=yL1a(xL-xLimit1+1);end
              if yL2a(xL-xLimit1+1)<yLimitTempL2;yLimitTempL2=yL2a(xL-xLimit1+1);end
              if yLimitTempR1>yLimitTempL1;yT=yLimitTempR1;else yT=yLimitTempL1;end
              if yLimitTempR2>yLimitTempL2; yB=yLimitTempL2;else yB=yLimitTempR2;end 
              if (xR-xLimit1+1)*(yB-yT+1)<LR(1)
                break;
              end              
              area=(xR-xL+1)*(yB-yT+1);
              if area>AreaNow
                AreaNow=area;
                if area>LR(1)
                  LR(1:5)=[area,xL,xR,yT,yB];
                end
              end             
            end%for xL=xNow:-1:xLimit1
          end%for xR=xNow:xLimit2          
        end%if xNow==xLimit1   
      end%if xLimit2-xLimit1>yLimit2-yLimit1
    end%if widthMax*heightMax>LR(1)
  end%for i=1:size(border,1)
  LR(9)=AreaNow;
end
function plotResult(ImBW,LR,x,y)
  hFigPlotLRC = findobj( 'Type', 'Figure', 'Tag', 'Fig$PlotLRC' );
  if isempty(hFigPlotLRC)
    screensize=get(0, 'MonitorPositions');
    hFigPlotLRC=figure('Tag','Fig$PlotLRC','Name','Largest shape',...
      'OuterPosition',[825,40,screensize(1,3)-825,screensize(1,4)-60]);
  end
  figure(hFigPlotLRC);
  cla reset;
  imsize=size(ImBW);
  L2=LR(3)-LR(2)+1;
  L1=LR(5)-LR(4)+1;
  imshow(ImBW,'border', 'tight','InitialMagnification','fit');
  axis on;
  hold on;
  plot([x,x(1)],[y,y(1)],'r-','LineWidth',1);%plot LR
  u=zeros(1,4);
  v=zeros(1,4);
  %draw arrows for x ,y
  uvL=max(abs([x(1)-x(2),x(1)-x(3),y(1)-y(2),y(1)-y(3)]))/8;%quiver length
  if x(3)==x(1)
    v(1)=-uvL;
    v(3)=uvL;
  else
    m1=((y(3)-y(1))/(x(3)-x(1)));
    u(1)=2*uvL/sqrt(1+m1^2);%make it longer by 2
    if x(1)<x(3)
      u(1)=-u(1);
    end  
    u(3)=-u(1);
    v(1)=m1*u(1);
    v(3)=-v(1);
  end
  if x(4)==x(2)
    u(2)=uvL;
    u(4)=-uvL;
  else
    m2=((y(4)-y(2))/(x(4)-x(2)));
    u(2)=uvL/sqrt(1+m2^2);
    if x(2)<x(4)
      u(2)=-u(2);
    end  
    u(4)=-u(2);
    v(2)=m2*u(2);
    v(4)=-v(2);
  end
  xa=x-u;
  ya=y-v;
  quiver(xa,ya,u,v,0,'MaxHeadSize', .2,'color',[.8 .6 .6]);
  xa=xa-0.2*u;
  ya=ya-0.2*v;
   text(xa,ya,{[num2str(x(1), '%.1f') ', ' num2str(y(1),'%.1f')],...
     [num2str(x(2), '%.1f') ', ' num2str(y(2),'%.1f')],...
     [num2str(x(3), '%.1f') ', ' num2str(y(3),'%.1f')],...
     [num2str(x(4), '%.1f') ', ' num2str(y(4),'%.1f')]},...
     'Color','r','HorizontalAlignment','center');
  text(x(1)+(x(3)-x(1))/2,y(1)+(y(3)-y(1))/2,...
    sprintf(['Area: ' num2str(L1) 'x'  num2str(L2) 'px=' num2str(LR(1)) 'px \n'...
    num2str(LR(6),'%.2f') ' degrees']),'HorizontalAlignment','center');
  % fprintf('Area: %.0fx%.0f px=%.0f px at %.2f degrees  \ncorner locations x, y: %.1f, %.1f | %.1f, %.1f | %.1f, %.1f | %.1f, %.1f\n',...
  % L1,L2,LR(1),LR(6),x(1),y(1),x(2),y(2),x(3),y(3),x(4),y(4));%Output to command window
  %angle
  phi=LR(6)/180*pi;
  L43=sqrt((x(3)-x(4))^2+(y(4)-y(3))^2);
  L32=sqrt((x(2)-x(3))^2+(y(3)-y(2))^2);
  if L43>L32
    %left is long side of LR
    dxT=(x(3)-x(4))/2;
    dyT=(y(3)-y(4))/2;
    txt=sprintf('%.2f °',LR(6));
    r=-sqrt(dxT^2+dyT^2);
    if r<-x(3)+imsize(2)*.1+1
      r=-x(3)+imsize(2)*.01+0.5;
      ytxt=y(3)+r*sin(phi/2);%+imsize(1)*.03;
      xtxt=x(3)+r*cos(phi/2)+imsize(2)*.01;
      text(xtxt,ytxt,txt,'Color','green','HorizontalAlignment','left','VerticalAlignment','top'); 
    else
      ytxt=y(3)-dyT/2;%txt y pos
      xtxt=x(3)+r*cos(phi/2)-imsize(2)*.01;
      text(xtxt,ytxt,txt,'Color','green','HorizontalAlignment','right','VerticalAlignment','baseline');
    end
    phic=0:0.02*phi:phi;
  else
    %right
    dxT=(x(2)-x(3))/2;
    dyT=(y(3)-y(2))/2;
    txt=sprintf('%.2f °',90-LR(6));
    r=sqrt(dxT^2+dyT^2);
    if r>imsize(2)*.9-x(3)-1;
      r=imsize(2)*.99-x(3);
      ytxt=y(3)-r*sin((pi/2-phi)/2); 
      xtxt=x(3)+r*cos((pi/2-phi)/2);
      text(xtxt,ytxt,txt,'Color','green','HorizontalAlignment','right','VerticalAlignment','top');    
    else
      ytxt=y(3)-dyT/2;
      xtxt=x(3)+r*cos((pi/2-phi)/2)+imsize(2)*.01;
      text(xtxt,ytxt,txt,'Color','green','HorizontalAlignment','left','VerticalAlignment','baseline');
    end  
    phic=-(0:0.02*(pi/2-phi):pi/2-phi);
   end
  plot([x(3)+r,x(3)],[y(3),y(3)],'g:','LineWidth',1);
  xc=x(3)+r*cos(phic);
  yc=y(3)+r*sin(phic);  
  plot(xc,yc,'g-','LineWidth',1);
end
end
