function [colors ] = colorGen(colorBegin, colorEnd, numColor)
%generate range of color varied from colorBegin to colorEnd

colors = [linspace(colorBegin(1),colorEnd(1),numColor)', linspace(colorBegin(2),colorEnd(2),numColor)', linspace(colorBegin(3),colorEnd(3),numColor)'];


end

