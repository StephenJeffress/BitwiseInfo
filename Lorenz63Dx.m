function [ dX ] = Lorenz63Dx(X,a)

%a(1) = sigma
%a(2) = rho
%a(3) = Beta
dX = zeros(3,1);
dX(1) = a(1)*(X(2)-X(1));
dX(2) = X(1)*(a(2)-X(3)) - X(2);
dX(3) = X(1)*X(2) - a(3)*X(3);

end

