function [ coefSol ] = nullMDS(m)
%computes the nullspace of matrix m using the Minimal Determining Set
%algorithm (Algorithm 1 in DOI: 10.1016/j.camwa.2016.10.014)

I_tilde = [];
v_kd = size(m,2);
I_bar = 1:size(m,2);
tol_sol = 1e-8; %tolerance for residual for solution

for i=1:v_kd
    %disp('enter loop')
    b_j = zeros(v_kd,1);
    b_j(i) = 1;
    red_rhs = m*b_j;
    setColumns = union(I_bar, I_tilde);
    unSetColumns = setdiff(1:v_kd,setColumns);
    %check if the constrained system has a solution
    tempSol = m(:,unSetColumns)\-red_rhs;
    residual = m(:,unSetColumns)*tempSol+red_rhs;
    if norm(residual)<tol_sol
        %disp('enter sol')
        b_j(unSetColumns)=tempSol;
        [~,max_br] = max(abs(b_j));
        I_tilde=union(I_tilde,max_br(1))  ;
    end
    I_bar=setdiff(I_bar,i);
end
minimal_determining_set=I_tilde;

%set of indices to solve for 
solIndex=setdiff(1:v_kd,I_tilde);
coefSol=zeros(v_kd,length(I_tilde));

for i=1:length(I_tilde)
    lhs=m;
    lhs(:,I_tilde)=[];
    rhs=-m(:,I_tilde(i));
    sol=lhs\rhs;
    tempB=zeros(1,v_kd);
    tempB(I_tilde(i))=1;
    tempB(solIndex)=sol;
    coefSol(:,i)=tempB;
end

end

