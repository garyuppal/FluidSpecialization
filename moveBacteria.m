   

function bact = moveBacteria(xmin,xmax,xbnd,ymin,ymax,ybnd,dt,...
        db,bact,veltype,vmax,rotation,vrad,delay,sticky)
    
    Nb = length(bact(bact(:,3)>=0)); % number alive 
    % (should be sorted to keep alive at top)
    % bact = (bi:)[xpos ypos s1]
    
    % === GET VELOCITY ===
    if veltype <=0 || sticky == 1 || delay == 1 
        bvel = 0;
    elseif veltype == 1
        maxrad = max(abs(ymin),ymax);
        bvx = (vmax/maxrad)*bact(1:Nb,2);
        bvy = zeros(Nb,1);
        bvel = [bvx bvy];
    elseif veltype == 2
        maxrad = max(abs(ymin),abs(ymax));
        bvx = vmax*(1 - ((bact(1:Nb,2)).^2)/(maxrad^2) );
        bvy = zeros(Nb,1);
        bvel = [bvx bvy];
    elseif veltype == 3
        [bvx, bvy] = brankine(rotation,vrad,bact(1:Nb,1),bact(1:Nb,2));

        bvel = [bvx, bvy];
    elseif veltype == 4
        bvx = vmax*ones(Nb,1);
        bvy = zeros(Nb,1);
        bvel = [bvx, bvy];
    else
        bvel = 0;
    end
    
    % === GET DIFFUSION ===
    stepsize = sqrt(db*4*dt);

    dpos = (rand(Nb,2)<0.5);
    dpos = stepsize*(dpos + (dpos(:,1)&dpos(:,2))*([-2 -1]) + ...
        (~(dpos(:,1)|dpos(:,2)))*([0 -1]));
    
    % === MOVE BACTERIA ===
    bact(1:Nb,1:2) = bact(1:Nb,1:2) + dpos + bvel*dt;
    
    % === CHECK BOUNDARY CONDITIONS ===
    tooSmallX = bact(1:Nb,1) < xmin;
    tooLargeX = bact(1:Nb,1) > xmax;    
    tooSmallY = bact(1:Nb,2) < ymin;
    tooLargeY = bact(1:Nb,2) > ymax;

    if xbnd == 0
        bact(1:Nb,1) = bact(1:Nb,1) + tooSmallX.*(xmax - xmin);
        bact(1:Nb,1) = bact(1:Nb,1) - tooLargeX.*(xmax - xmin);
    elseif xbnd == 1
        bact(1:Nb,1) = bact(1:Nb,1) + 2*tooSmallX.*(xmin-bact(1:Nb,1));
        bact(1:Nb,1) = bact(1:Nb,1) - 2*tooLargeX.*(bact(1:Nb,1)-xmax);
    elseif xbnd == 2
        bact(1:Nb,3) = bact(1:Nb,3) - tooLargeX.*(100 + bact(1:Nb,3));
        bact(1:Nb,3) = bact(1:Nb,3) - tooSmallX.*(100 + bact(1:Nb,3)); % kill bacteria that escape
    end
    if ybnd == 0
        bact(1:Nb,2) = bact(1:Nb,2) + tooSmallY.*(ymax - ymin);
        bact(1:Nb,2) = bact(1:Nb,2) - tooLargeY.*(ymax - ymin);
    elseif ybnd == 1
        bact(1:Nb,2) = bact(1:Nb,2) + 2*tooSmallY.*(ymin-bact(1:Nb,2));
        bact(1:Nb,2) = bact(1:Nb,2) - 2*tooLargeY.*(bact(1:Nb,2)-ymax);
    elseif ybnd == 2
        bact(1:Nb,3) = bact(1:Nb,3) - tooLargeY.*(100 + bact(1:Nb,3));
        bact(1:Nb,3) = bact(1:Nb,3) - tooSmallY.*(100 + bact(1:Nb,3)); % kill bacteria that escape
    end

end