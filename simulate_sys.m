function xy = simulate_sys( x, y, u, delta_t)

% Return [x1 x2 x3 y]

% Evaluate the derivatives
x_dot(1) = -x(1)+u*(2+x(3)^2)/(1+x(3)^2);
x_dot(2) = x(3);
x_dot(3) = x(1)*x(3)+u;

y_dot = x_dot(2);

% Integrate to get a new x,y
x(1) = x(1)+x_dot(1)*delta_t;
x(2) = x(2)+x_dot(2)*delta_t;
x(3) = x(3)+x_dot(3)*delta_t;
y = y+y_dot*delta_t;

xy = [x y];
