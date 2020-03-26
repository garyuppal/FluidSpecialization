function F = displayGraphics(xmin,xmax,Nx,ymin,ymax,Ny,bact,alive,...
            c1,c2,cw,t,figHand)
        
    % colors:
    genc = (1/300)*[169 209 141];
    yelc = (1/255)*[255 208 59];
    redc = (1/255)*[197 90 17]; 
    chtc = (1/255)*[56 87 35];

    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;
    x = xmin:dx:xmax; % including ghost nodes x(1) and x(Nx+3)
    y = ymin:dy:ymax; % including ghost nodes y(1) and y(Ny+3)
    
    
    if alive == 1 
        Nb = length(bact(bact(:,3)>=0)); % number alive 
        colors = zeros(Nb,3);
        xpos = bact(1:Nb,1);
        ypos = bact(1:Nb,2);
        for bi = 1:Nb
            colors(bi,:) = genc.*( bact(bi,3) > 0).*( bact(bi,4) > 0) ...
                + redc.*( bact(bi,3) > 0).*( bact(bi,4) == 0) ...
                + yelc.*( bact(bi,3) == 0).*( bact(bi,4) > 0) ...
                + chtc.*( bact(bi,3) == 0).*( bact(bi,4) == 0);
            if bact(bi,3) <0
                disp('problem!...');
                disp(bact(bi,3));
            end
        end
    end
    
    mc1 = 1.0/(max(max(c1)) + 0.0001);
    mc2 = 1.0/(max(max(c2)) + 0.0001);
    mcw = 1.0/(max(max(cw)) + 0.0001);
    clrscl = max([mc1 mc2 mcw]); % or mean...
    
    Clrs = zeros(Ny+1,Nx+1,3);
    for xi = 2:Nx+2
        for yi = 2:Ny+2
            Clrs(yi-1,xi-1,:) = 2*mc1*redc.*c1(xi,yi) ...
                + 2*mc2*yelc.*c2(xi,yi) + [0 0 clrscl]*cw(xi,yi);
        end
    end        
    
    figure(figHand);
    image(x,y,Clrs);
    set(gca,'YDir','normal');
    shading interp;
    hold on;
    if alive ==1
        scatter(xpos,ypos,10,colors,'filled'); %,'^',...
%             'MarkerEdgeColor',(1/256)*[192 192 192],'LineWidth',0.5);
    end
    axis([xmin xmax ymin ymax]);

    title(sprintf('time = %1.8f',t)); 
    pause(0.0001);
    hold off;
    
    F = getframe(gcf);

end

