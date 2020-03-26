% secrete onto chemical grid:

function [cp1, cp2, cpw] = bsecrete(xmin,xmax,Nx,dx,ymin,ymax,Ny,dy,bact,sw)

ctmp1 = zeros(Nx+3,Ny+3);
ctmp2 = zeros(Nx+3,Ny+3); % including ghost node
ctmpw = zeros(Nx+3,Ny+3); % including ghost node

invArea = (dx*dy); 

Nb = length(bact(bact(:,3)>=0)); % number alive 

for bi = 1:Nb
    xp = bact(bi,1);
    yp = bact(bi,2);
    s1 = bact(bi,3);
    s2 = bact(bi,4);
    
    [xi,yi] = findGridIndex(xmin,xmax,Nx,ymin,ymax,Ny,xp,yp);
    ctmp1(xi,yi) = ctmp1(xi,yi) + s1*invArea;
    ctmp2(xi,yi) = ctmp2(xi,yi) + s2*invArea;
    ctmpw(xi,yi) = ctmpw(xi,yi) + sw*invArea;
end

cp1 = ctmp1;
cp2 = ctmp2;
cpw = ctmpw;
end