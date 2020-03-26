function [data, dataChem] = simulator(xmin,xmax,xbnd,Nx,ymin,ymax,ybnd,Ny,tmax,dt,...
    xbmin,xbmax,ybmin,ybmax, ...
    db,d1,d2,dw,a1,a2,aor,a12,aw,b1,b2,k1,k2,kor,k12,kw,s1,s2,sw,l1,l2,lw, ...
    gamma1,gamma2,delta, ...
    s1rand,s2rand,sbinary, ...
    veltype,vmax,rotation,vrad,...
    mutProb1,mutDiff1,mdiffrand1,mutProb2,mutDiff2,mdiffrand2,binaryMutation,...
    sticky,initialNumBacteria,initialNumGroups,aligngroups,maxBacteria,...
    velDelay,mutDelay,repDelay,...
    savePeriod,saveChemicals,...
    graphics,framePeriod,saveVid,vidOption,vidFile,...
    tsave,dataFile)
%------------------------------------------------------------------------

% ===== INITIALIZE VARIABLES =====
rng('shuffle'); % seed random number generator
alive = 1;      % set state to living
ogs1 = s1;      % for binary mutation: initial secretion rate
ogs2 = s2;

% === DISCRETIZE SPACE === 
dx = (xmax - xmin)/Nx;
dy = (ymax - ymin)/Ny;
x = (xmin-dx):dx:(xmax+dx); % including ghost nodes x(1) and x(Nx+3)
y = (ymin-dx):dy:(ymax+dx); % including ghost nodes y(1) and y(Ny+3)

Nt = ceil(tmax/dt);         % number of timesteps
t = 0.0;                    % initialize time
tsaveInd = ceil(tsave/dt);

% === INITILIZE SAVE VARIABLE ===
nsave = ceil(Nt/savePeriod);
dataRec = cell(1,nsave);
if saveChemicals == 1
    dataC = cell(1,nsave);
end % if also saving chemicals
svi = 0;

% === GENERATE VELOCITY FIELD ===
[velX,velY] = genVelocity(veltype,vmax,rotation,vrad,x,y);

% === INITIALIZE BACTERIA AND CHEMICAL FIELDS ===
bact = constructBacteria(xbmin,xbmax,ybmin,ybmax,s1,s2,initialNumBacteria,...
    initialNumGroups,aligngroups,maxBacteria,s1rand,s2rand,sbinary); 

% initialize chemicals:
c1 = zeros(Nx+3,Ny+3);
c2 = zeros(Nx+3,Ny+3); % including ghost nodes
cw = zeros(Nx+3,Ny+3); % including ghost nodes

% === GRAPHICS ===
if saveVid == 1 && graphics == 1
    if vidOption == 2
        vid = VideoWriter(vidFile,'Uncompressed AVI');
    else
        vid = VideoWriter(vidFile,'MPEG-4');
    end
    open(vid);
end

% figure for video
if graphics == 1
    scrsz = get(groot,'ScreenSize');
    % scrsz = [1 1 width height]
    % position: [left bottom width height]
    h = figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]);
end


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% ===== LOOP THROUGH TIME =====
for ti = 1:Nt
    % === BACTERIA DYNAMICS ===
    bact = moveBacteria(xmin,xmax,xbnd,ymin,ymax,ybnd,dt,...
        db,bact,veltype,vmax,rotation,vrad,(t<velDelay),sticky);
    [cp1, cp2, cpw] = bsecrete(xmin,xmax,Nx,dx,ymin,ymax,Ny,dy,bact,sw);

    if t > repDelay
        [bact,alive] = reproduceBacteria(xmin,xmax,Nx,...
            ymin,ymax,Ny,dt,...
             a1,a2,aor,a12,aw,b1,b2,k1,k2,kor,k12,kw,...
             bact,c1,c2,cw,gamma1,gamma2,delta);
    end
    
    % possibly mutate bacteria:
    if t > mutDelay
        bact = mutateBacteria2(bact,...
             mutProb1,mutDiff1,mdiffrand1,...
             mutProb2,mutDiff2,mdiffrand2,binaryMutation,ogs1,ogs2);
    end
    
    if alive == 0
        disp('Everybody died!');
        break;
    end
    
    % === CHEMICAL DYNAMICS ===
    % apply boundary conditions:
    [c1, c2, cw] = updateGhosts(c1,c2,cw,Nx,Ny,xbnd,ybnd);
    
    % chem 1:
    cppdiff1 = d1*laplacian(dx,dy,c1);
    if veltype <= 0 || t < velDelay
        cppadv1 = 0;
    else
        cppadv1 = chemAdvect(dx,dy,velX,velY,c1);
    end

    c1 = c1 + (cp1 + cppdiff1 + cppadv1 - l1*c1)*dt;
    
    % chem 2:
    cppdiff2 = d2*laplacian(dx,dy,c2);
    if veltype <= 0 || t < velDelay
        cppadv2 = 0;
    else    
        cppadv2 = chemAdvect(dx,dy,velX,velY,c2);
    end
  
    c2 = c2 + (cp2 + cppdiff2 + cppadv2 - l2*c2)*dt;
    
    % chem w:
    cppdiffw = dw*laplacian(dx,dy,cw);
    if veltype <= 0 ||  t < velDelay
        cppadvw = 0;
    else    
        cppadvw = chemAdvect(dx,dy,velX,velY,cw);
    end
  
    cw = cw + (cpw + cppdiffw + cppadvw - lw*cw )*dt;
    
    % update time:
    t = t + dt;
    
    % === SAVE INTERMEDIATE ===
    if ismember(ti,tsaveInd)
        save_file = sprintf('T%.3f_%s',t,dataFile);
        parsave(save_file,dataRec);
    end % if saving

    
    % == SAVE DATA ==
    if mod(ti-1,savePeriod) == 0
        svi = svi + 1;
        Nb = length(bact(bact(:,3)>=0)); % number alive 
        dataRec{svi} = bact(1:Nb,:);
        if saveChemicals == 1
            dataC{svi} = [c1; c2; cw];
        end % if also saving chemicals
    end
    
    % == DISPLAY GRAPHICS ==
    if graphics == 1 && mod(ti-1,framePeriod) == 0
        F = displayGraphics(xmin,xmax,Nx,ymin,ymax,Ny,bact,alive,...
            c1,c2,cw,t,h);
        if saveVid == 1
            writeVideo(vid,F);
        end
    end

end % end loop over time
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

% ===== GENERATE OUTPUT =====
% close video:
if saveVid == 1 && graphics == 1
    close(vid);
end

% return data:
data = dataRec;
if saveChemicals == 1
    dataChem = dataC;
else
    dataChem = 0;
end

end % END OF FUNCTION

