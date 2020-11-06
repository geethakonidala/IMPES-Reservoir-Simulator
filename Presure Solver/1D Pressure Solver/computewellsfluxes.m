function q = computewellsfluxes(pw, PI, Lamda, cell, p)

    N = length(p);

    q = zeros(N,1);

for w = 1:length(PI)
    
    q(cell(w)) = q(cell(w),1)+PI(w)*Lamda(cell(w))*(pw(w) - p(cell(w)));   % check this line
end
end