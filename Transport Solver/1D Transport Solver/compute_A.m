function [A] = compute_A(N, T)
A = zeros(N, N);
    for i = 1:N
            if (i>1)            % Left interface exhist
         %T(i)*(p(i)-p(i-1))
                A(i,i)= A(i,i)+ T(i);
                A(i,i-1)= - T(i);
            end
            if (i<N)             % Right Interface exhist
          %T(i+1)*(p(i)- p(i+1))
                A(i,i)=  A(i,i)+ T(i+1); % Adding previous A(i,i) to this as well
                A(i,i+1)= - T(i+1);
            end
    end
end 