% function im = drawEllipse(area,epsilon,rot,fgcol,bgcol)
%
% Returns a matrix with an ellipse image, with specified size and
% eccenricity.
%
% INPUT
%  area      : ellipse size (pixels^2)
%  epsilon   : ellipse eccentricity
%  rot       : rotation (deg)
%  fgcol     : foreground color
%  bgcol     : background color
%
% OUTPUT
%  im        : image with the requested image
%
% Written by RvdB, Apr 2010

function im = drawEllipse(area,epsilon,rot,fgcol,bgcol)

rot = pi*(rot-90)/180;  

% compute lengths of long and short axis
b = sqrt(area * sqrt(1 - epsilon^2) / pi);
a = area / (pi*b);

% convert to diameters
d1 = round(2*b);
d2 = round(2*a);

% make sure that d1 is the minor axis
if (d1>d2)
    d3=d1;
    d1=d2;
    d2=d3;
end

% draw ellipse
im = ones(2*d2,2*d2)*bgcol;
minX = -d2;
maxX = minX + 2*d2 - 1; 
[X Y] = meshgrid(minX:maxX,minX:maxX);
X_new = X * cos(rot) - Y * sin(rot);
Y = X * sin(rot) + Y * cos(rot);
X = X_new;
idx = (X.^2/(d1/2)^2 + Y.^2/(d2/2)^2)<1;
im(idx) = fgcol;

% crop
while im(:,1)==bgcol
    im = im(:,2:end);
end
while im(1,:)==bgcol
    im = im(2:end,:);
end
while im(end,:)==bgcol
    im = im(1:end-1,:);
end
while im(:,end)==bgcol
    im = im(:,1:end-1);
end


