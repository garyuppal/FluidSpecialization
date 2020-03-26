

function  bact = constructBacteria(xmin,xmax,ymin,ymax,s1,s2,initialNumBacteria,...
    initialNumGroups,aligngroups,maxBacteria,s1rand,s2rand,sbinary)

% make bacteria into a matrix with columns: 
%[ xpos ypos s1 s2] -- can add more...

btemp = zeros(maxBacteria,4);

if initialNumGroups > 0
    if aligngroups == 1
        grpPosY = 0.5*(ymax-ymin)*ones(initialNumGroups,1) + ymin;
        if initialNumGroups == 1
            grpPosX = 0.5*(xmax-xmin)+xmin; % center singleton group
        else
            tmpspacing = linspace(xmin,xmax,initialNumGroups+1);
            grpPosX = transpose(tmpspacing(1:initialNumGroups));
        end
    else
        grpPosX = (xmax-xmin)*rand(initialNumGroups,1) + xmin;
        grpPosY = (ymax-ymin)*rand(initialNumGroups,1) + ymin;
    end % if aligning groups

    grpPos = [grpPosX grpPosY];
end % create group positions if needed

btemp(:,3) = -100;
for bi = 1:initialNumBacteria
    if initialNumGroups > 0
        j = mod(bi,initialNumGroups) + 1;
        p1 = grpPos(j,1);
        p2 = grpPos(j,2);
    else
        p1 = (xmax-xmin)*rand + xmin;
        p2 = (ymax-ymin)*rand + ymin;
    end
    
    % random and binary:
    if s1rand == 1
        s1act = s1*rand;
    else
        s1act = s1;
    end
    if s2rand == 1
        s2act = s2*rand;
    else
        s2act = s2;
    end
    if sbinary == 1
        pselect = rand;
        if pselect > 0.5
            s1act = 0;
        else
            s2act = 0;
        end
    end
    
    btemp(bi,1) = p1;           % initial x position
    btemp(bi,2) = p2;           % initial y position
    btemp(bi,3) = s1act;        % initial s1 - public good secretion
    btemp(bi,4) = s2act;        % initial s2 - public good secretion
end % for each bacteria

bact = btemp;
end


