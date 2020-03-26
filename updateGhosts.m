function   [c1, c2, cw] = updateGhosts(c1,c2,cw,Nx,Ny,xbnd,ybnd)

    % x boundary conditions:
    if xbnd == 0
        % periodic:
        c1(1,:) = c1(Nx+2,:); 
        c1(Nx+3,:) = c1(2,:); 
        c2(1,:) = c2(Nx+2,:); 
        c2(Nx+3,:) = c2(2,:); 
        cw(1,:) = cw(Nx+2,:); 
        cw(Nx+3,:) = cw(2,:); 
    elseif xbnd == 1
        % no flux:
        c1(1,:) = c1(2,:); 
        c1(Nx+3,:) = c1(Nx+2,:); 
        c2(1,:) = c2(2,:); 
        c2(Nx+3,:) = c2(Nx+2,:); 
        cw(1,:) = cw(2,:); 
        cw(Nx+3,:) = cw(Nx+2,:); 
    elseif xbnd == 2
        % no flux:
        c1(1,:) = 0; 
        c1(Nx+3,:) = 0; 
        c2(1,:) = 0; 
        c2(Nx+3,:) = 0; 
        cw(1,:) = 0; 
        cw(Nx+3,:) = 0; 
    else
        % periodic (default):
        c1(1,:) = c1(Nx+2,:); 
        c1(Nx+3,:) = c1(2,:); 
        c2(1,:) = c2(Nx+2,:); 
        c2(Nx+3,:) = c2(2,:); 
        cw(1,:) = cw(Nx+2,:); 
        cw(Nx+3,:) = cw(2,:); 
    end % x boundary conditions
    
    % y boundary conditions:
    if ybnd == 0
        % periodic:
        c1(:,1) = c1(:,Ny+2); 
        c1(:,Ny+3) = c1(:,2); 
        c2(:,1) = c2(:,Ny+2); 
        c2(:,Ny+3) = c2(:,2); 
        cw(:,1) = cw(:,Ny+2); 
        cw(:,Ny+3) = cw(:,2); 
    elseif ybnd == 1
        % no flux:
        c1(:,1) = c1(:,2); 
        c1(:,Ny+3) = c1(:,Ny+2); 
        c2(:,1) = c2(:,2); 
        c2(:,Ny+3) = c2(:,Ny+2);
        cw(:,1) = cw(:,2); 
        cw(:,Ny+3) = cw(:,Ny+2); 
   elseif ybnd == 2
        % no flux:
        c1(:,1) = 0; 
        c1(:,Ny+3) = 0; 
        c2(:,1) = 0; 
        c2(:,Ny+3) = 0;
        cw(:,1) = 0; 
        cw(:,Ny+3) = 0; 
    else 
        % no flux (default)
        c1(:,1) = c1(:,2); 
        c1(:,Ny+3) = c1(:,Ny+2); 
        c2(:,1) = c2(:,2); 
        c2(:,Ny+3) = c2(:,Ny+2); 
        cw(:,1) = cw(:,2); 
        cw(:,Ny+3) = cw(:,Ny+2); 
    end % y boundary conditions
end

