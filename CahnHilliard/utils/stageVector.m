function [C_alphaF,Cdot_alphaM,Cdot1] = stageVector( C0,Cdot0,C1,alphaM,alphaF,gamma,timeStep )
%stage Vector CahnHilliard (eq 27.1 & 27.2);

[Cdot1] = VecWAXPY(-1.0,C0,C1 );
[Cdot1] = VecAXPBY(Cdot1,1-1/gamma,1/(gamma*timeStep),Cdot0 );

[C_alphaF] = VecWAXPY(-1.0,C0,C1 );
[C_alphaF] = VecAYPX( C_alphaF,alphaF,C0);

[Cdot_alphaM] = VecWAXPY(-1.0,Cdot0,Cdot1 );
[Cdot_alphaM] = VecAYPX( Cdot_alphaM,alphaM,Cdot0);

end
