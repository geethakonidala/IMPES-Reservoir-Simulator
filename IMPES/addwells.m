function [A,q] = addwells(A,q,lamda_total,pw,PI,cell)
% intialization of parameters
for i=1:length(pw) %d number of wells
    A(cell(i),cell(i))= A(cell(i),cell(i))+(PI(i)*lamda_total(cell(i)));
    q(cell(i))        = q(cell(i))+ PI(i)*lamda_total(cell(i))*pw(i);
end

end 