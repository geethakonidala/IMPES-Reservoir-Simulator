function [q ] = ComputeWellsFluxes( P,PI,lambda,cells,pw,N)
q=zeros(N,1);
n_wells=length(pw);
for w=1:n_wells
    q(cells(w))=PI(w)*lambda(cells(w))*(pw(w)-P(cells(w)));
end
end

