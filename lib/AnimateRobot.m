function AnimateRobot(H,AAA,P,h,sd,plotpath)
    % H - graphic handles of objects to draw
    % AAA - obtained from CalculateRobotMotion
    % P - object to represent
    % h - robot's links' graphic handle
    % sd - pause between frames, in seconds
    % plotpath - boolean, draw point on each position of hand over time
    
    if nargin < 6
        plotpath = false;
    end
    
    M = size(AAA, 4);  % number of positions -> nº of frames
    
    for m = 1:M
        
        Org = LinkOrigins(AAA(:,:,:,m));
        h.XData = Org(1,:);
        h.YData = Org(2,:);
        h.ZData = Org(3,:);
    
        Ak = eye(4);
        for n = 1:size(AAA, 3)
            Ak = Ak*AAA(:,:,n,m);
            Pn = Ak*P;
            H{n}.Vertices = Pn(1:3,:)';
        end
    
        if plotpath
            plot3(Ak(1,4), Ak(2,4), Ak(3,4), '.r');
        end
    
        pause(sd)
    end

end

