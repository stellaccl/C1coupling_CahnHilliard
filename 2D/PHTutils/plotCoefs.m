function   plotCoefs( coefs ,side )
%plot coefs for surfaces
switch side
    case 1
        temCoefs=coefs(:,:,1);
        
        for i=1:size(temCoefs,2)
            p=temCoefs(:,i);
            hold on
            plot3(p(1),p(2),p(3),'*k')
        end
        
    case 2
        temCoefs=coefs(:,end,:);
        
        for i=1:size(temCoefs,3)
            p=temCoefs(:,1,i);
            hold on
            plot3(p(1),p(2),p(3),'*k')
            
        end
    case 3
        temCoefs=coefs(:,:,end);
        
        for i=1:size(temCoefs,2)
            p=temCoefs(:,i);
            hold on
            plot3(p(1),p(2),p(3),'*k')
        end
    case 4
        temCoefs=coefs(:,1,:);
        
        for i=1:size(temCoefs,3)
            
            p=temCoefs(:,1,i);
            hold on
            plot3(p(1),p(2),p(3),'*k')
        end
end


end

