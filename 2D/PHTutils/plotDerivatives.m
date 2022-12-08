function  [up2,vp2]=plotDerivatives(nurbsInfo,edge,npts,color )
%plot the x and y derivaitves on surfaces

tt = linspace(0.0,1.0,npts);

switch edge
    case 1
        plotLoc={tt,0 };
    case 2
        plotLoc={1,tt };
    case 3
        plotLoc={tt,1 };
    case 4
        plotLoc={0,tt };  
end


dsrf = nrbderiv (nurbsInfo);

[p1, dp] = nrbdeval(nurbsInfo, dsrf, plotLoc);

up2 = vecnorm(dp{1})
vp2 = vecnorm(dp{2})

hold on;
plot3(p1(1,:),p1(2,:),p1(3,:),'ro');
h1 = quiver3(p1(1,:),p1(2,:),p1(3,:),up2(1,:),up2(2,:),up2(3,:));
%h2 = quiver3(p1(1,:),p1(2,:),p1(3,:),vp2(1,:),vp2(2,:),vp2(3,:));
set(h1,'Color',color);
%set(h2,'Color',color);


end

