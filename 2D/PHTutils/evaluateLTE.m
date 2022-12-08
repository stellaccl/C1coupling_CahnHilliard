function [error] = evaluateLTE( CBE,CLTE,C0,Cnew,prevSol,step,timeStep,time,previousTime )
if step==0
    
    X=Cnew;
    Y=X+CLTE;
    diff=(X-Y);
    error=norm(diff);
else
    X=Cnew;
    h=timeStep;
    h_prev=time-previousTime;
    
    a=1+h_prev/h;

    scal=[1/a;-1/(a-1);1/(a*(a-1))];
    vecs=[Cnew,C0,prevSol];
    tempX=zeros(size(Cnew));
    
    for i=1:3
        tempX=tempX+scal(i)*vecs(:,i);
    end

    Y=X;
    Y=Y+tempX;
    diff=(X-Y);
    error=norm(diff);
end

end

