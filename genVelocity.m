% generate velocity function:

function [vx,vy] = genVelocity(veltype,vmax,rotation,vrad,x,y)
% initialize velocity:
% veltype options:
% 1 = couette
% 2 = poisuille
% 3 = rankine
% 4 = constant 

Nx = length(x);
Ny = length(y);

vxtmp = zeros(Nx,Ny);
vytmp = zeros(Nx,Ny);

PI = 3.1415926535897931;

if veltype == 1
    for yi = 1:Ny
        maxrad = max(y(end-1),y(2));
        vxtmp(:,yi) = (vmax/maxrad)*y(yi)*ones(Nx,1);
    end
elseif veltype == 2
    for yi = 1:Ny
        maxrad = max(abs(y(end-1)),abs(y(2)));
        vxtmp(:,yi) = vmax*( 1 - ((y(yi))^2)/(maxrad^2) )*ones(Nx,1);
    end
elseif veltype == 3
    for i = 1:Ny
        for j = 1:Nx
            theta = atan2(y(i),x(j));
            dist = sqrt(y(i)^2 + x(j)^2);
            if dist <= vrad
                inside = (rotation*dist)/(2*PI*vrad^2);
                vxtmp(j,i) = (sin(theta))*inside;
                vytmp(j,i) = -cos(theta)*inside;
            else
                outside = rotation/(2*PI*dist);
                vxtmp(j,i) = (sin(theta))*outside;
                vytmp(j,i) = -cos(theta)*outside;
            end
        end
    end
elseif veltype == 4
    vxtmp(:,:) = vmax;
end

vx = vxtmp;
vy = vytmp;
end