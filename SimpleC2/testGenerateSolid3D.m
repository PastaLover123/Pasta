params = createParams();

x = (0:(params.L-1))*params.dx;
[X,Y,Z]=meshgrid(x,x,x);

q1 = 1/sqrt(2)*[1 1 0];
q2 = 1/sqrt(2)*[1 0 1];
q3 = 1/sqrt(2)*[0 1 1];
q4 = q1-q2;
q5 = q2-q3;
q6 = q3-q1;

% solid = cos(X).*cos(Y)+cos(X).*cos(Z)+cos(Z).*cos(Y);
solid = cos(q1(1)*X+q1(2)*Y+q1(3)*Z) + cos(q2(1)*X+q2(2)*Y+q2(3)*Z) + cos(q3(1)*X+q3(2)*Y+q3(3)*Z) +...
        cos(q4(1)*X+q4(2)*Y+q4(3)*Z) + cos(q5(1)*X+q5(2)*Y+q5(3)*Z) + cos(q6(1)*X+q6(2)*Y+q6(3)*Z);

isosurface(X,Y,Z,solid,4)