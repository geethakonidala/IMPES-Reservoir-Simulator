function Index = computeIndex(i, j, Nx)
Index(I, J) = (j-1)*Nx +1

        I = 1+rem(Index - 1,Nx);
        J = floor((Index - 1)/Nx + 1);

end 