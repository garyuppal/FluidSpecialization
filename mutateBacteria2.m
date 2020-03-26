

function bact = mutateBacteria2(bact,mprob1,mdiff1,mdr1,...
            mprob2,mdiff2,mdr2,binaryMutation,ogs1,ogs2)
        
    Nb =  length(bact(bact(:,3)>=0)); % number alive 
    
    prob1 = rand(Nb,1);
    prob2 = rand(Nb,1);
    mutate1 = prob1<mprob1;
    mutate2 = prob2<mprob2;
    
    % mutation strengths:
    mu = 0;  
    if mdr1 == 1
        ms1 = mdiff1*(2*rand(Nb,1) - 1);
    else
        ms1 = normrnd(mu,mdiff1,[Nb 1]); % gaussian random walk
    end
    if mdr2 == 1
        ms2 = mdiff2*(2*rand(Nb,1) - 1);
    else
        ms2 = normrnd(mu,mdiff2,[Nb 1]); % gaussian random walk
    end
    
    if binaryMutation == 1
        bact(1:Nb,3) = bact(1:Nb,3) + ogs1*(mutate1.*( (bact(1:Nb,3)==0) - ...
            (bact(1:Nb,3)==ogs1) )); 
        bact(1:Nb,4) = bact(1:Nb,4) + ogs2*(mutate2.*( (bact(1:Nb,4)==0) - ...
            (bact(1:Nb,4)==ogs2) )); 
    else
        nds1 = bact(1:Nb,3) + mutate1.*ms1;
        nds2 = bact(1:Nb,4) + mutate2.*ms2;
        bact(1:Nb,3) = nds1.*(nds1>0); %abs(bact(1:Nb,3) + mutate1.*ms1);
        bact(1:Nb,4) = nds2.*(nds2>0); %abs(bact(1:Nb,4) + mutate2.*ms2); 
        % absorbing lower bc
    end

end
