   

function cpp = chemAdvect(dx,dy,velX,velY,c)

    Nx = size(c,1) - 3;
    Ny = size(c,2) - 3;
    mmx = (1/dx);
    mmy = (1/dy);
    
    % including ghosts there are N+3 nodes along each dimension
    ctmp = zeros(Nx+3,Ny+3);

for xi = 2:(Nx+2)
    for yi = 2:(Ny+2)
        if velX(xi,yi) >= 0 && velY(xi,yi) >= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi,yi)-c(xi-1,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi)-c(xi,yi-1))*mmy;
        elseif velX(xi,yi) >= 0 && velY(xi,yi) <= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi,yi)-c(xi-1,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi+1)-c(xi,yi))*mmy;
        elseif velX(xi,yi) <= 0 && velY(xi,yi) >= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi+1,yi)-c(xi,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi)-c(xi,yi-1))*mmy;
        elseif velX(xi,yi) <= 0 && velY(xi,yi) <= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi+1,yi)-c(xi,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi+1)-c(xi,yi))*mmy;
        end
    end
end

cpp = ctmp;
end

