function [refPoint] = assignRefPoint( ic, uref,vref,edgeA,edgeB )
%assign reference point based on edge

switch edgeA
    case 1
        urefA=uref(ic);
        vrefA=-1;
    case 2
        urefA=1;
        vrefA=vref(ic);
    case 3
        urefA=uref(ic);
        vrefA=1;
    case 4
        urefA=-1;
        vrefA=vref(ic);
end

switch edgeB
    case 1
        urefB=uref(ic);
        vrefB=-1;
    case 2
        urefB=1;
        vrefB=vref(ic);
    case 3
        urefB=uref(ic);
        vrefB=1;
    case 4
        urefB=-1;
        vrefB=vref(ic);
end

refPoint=struct;
refPoint.urefA=urefA;
refPoint.vrefA=vrefA;
refPoint.urefB=urefB;
refPoint.vrefB=vrefB;



end

