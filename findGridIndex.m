

function  [xi,yi] = findGridIndex(xmin,xmax,Nx,ymin,ymax,Ny,xp,yp)
dx = (xmax - xmin)/Nx;
dy = (ymax - ymin)/Ny;

x1 = xp - xmin + dx;
y1 = yp - ymin + dy; % add ghost node distance

x2 = abs(x1 * (1.0/dx) - 0.000000000001);
y2 = abs(y1 * (1.0/dy) - 0.000000000001);

xi = floor(x2) + 1;
yi = floor(y2) + 1;
end


