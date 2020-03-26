function  [bact,state] = reproduceBacteria(xmin,xmax,Nx,ymin,ymax,Ny,dt,...
        a1,a2,aor,a12,aw,b1,b2,k1,k2,kor,k12,kw,...
        bact,c1,c2,cw,gamma1,gamma2,delta)

    maxBacteria = size(bact,1);
    Nb =  length(bact(bact(:,3)>=0)); % number alive 
    
    offspring = zeros(1,Nb); % number of copies (or deaths) to produce

    % ===  DETERMINE NUMBER OF OFFSPRING ===
    for bi = 1:Nb
        % get fitness for each bacteria:
        xp = bact(bi,1);
        yp = bact(bi,2);
        s1 = bact(bi,3);
        s2 = bact(bi,4);

        [xi,yi] = findGridIndex(xmin,xmax,Nx,ymin,ymax,Ny,xp,yp);
        fit = findFitness(xi,yi,c1,c2,cw,a1,a2,aor,a12,aw,b1,b2,k1,k2,kor,k12,kw,...
            s1,s2,gamma1,gamma2,delta)*dt;
        
        fitInt = floor(fit);
        fitDec = fit - fitInt;

        % reproduce with probability prepr 
        prepr = rand;
        fitProb = 0;
        if prepr < fitDec
            fitProb = 1;
        end
        % assign number of offspring:
        if fitInt >= 0
            offspring(bi) = fitInt + fitProb;
        else
            offspring(bi) = -1;
        end 
    end % for each bacteria -- record number of offspring to produce

    totalRep = sum(offspring);
    living = sum(offspring>=0);

    if living <= 0
        state = 0;
        return;
    elseif totalRep >= maxBacteria
        disp('Too many bacteria!');
        state = 0;
        return;
    end %if bacteria go extinct or exceed maximum
    
    % === MAKE OFFSPRING ===
    % first ``kill'' to make room:
    bact(offspring<0,3) = -100; 
    % push dead to end:
    [~,ind] = sort(squeeze(bact(:,3)),'descend');
    bact = bact(ind,:);
    
    start = find(bact(:,3)<0,1); % start after last living bacteria
    
    tempind = 1:length(ind);
    ogind(ind) = tempind;
    
    % == REPRODUCE ==
    k = start;
    
    for bi = 1:Nb
        if offspring(bi) > 0
            for repj = 1:offspring(bi)
                bact(k,:) = bact(ogind(bi),:); % copy original position
                k = k + 1;
            end
        end % for offspring
    end % can posssibly speed this up by not looping -- matrix operation...
    
    
  state = 1;  
  end % END OF FUNCTION
  
  