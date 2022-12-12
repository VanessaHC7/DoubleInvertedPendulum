function N_Inverted_Pendulum_V3(N)

% Term project 2.151

% Modelling solver (Lagrange Equations)

syms t
syms q dq ddq [(N+1) 1]

q = q;
dq = dq;
ddq = ddq;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

syms m I l b [(N+1) 1] real
syms g f

mv = m;
Iv = I;
lv = l;
bv = b;

lv(1) = 0;
Iv(1) = 0;

% Rotation matrix function

R = @(theta) [cos(pi/2-theta) -sin(pi/2-theta); sin(pi/2-theta) cos(pi/2-theta)];

syms x y [(N+1) 1]

x(1) = q(1);
y(1) = 0;

for i = 2:(N+1)
    x(i) = x(i-1) + [1 0]*R(q(i-1))*[lv(i-1)/2 0]' + [1 0]*R(q(i))*[lv(i)/2 0]';
    y(i) = y(i-1) + [0 1]*R(q(i-1))*[lv(i-1)/2 0]' + [0 1]*R(q(i))*[lv(i)/2 0]';
end

clear T V L Q Eqn
syms Q Eqn [(N+1) 1]
syms L 

L = 0;
Q(:) = 0;
Q(1) = f;

for i = 1:(N+1)
    T(i) = (1/2)*mv(i)*ddt(x(i))^2 + (1/2)*mv(i)*ddt(y(i))^2 + (1/2)*Iv(i)*ddt(q(i))^2;
    V(i) = mv(i)*g*y(i);
    L = L + T(i) - V(i);
    Q(i) = Q(i) - bv(i)*dq(i);
end

Eqn = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q

Mn = simplify(jacobian(Eqn,ddq))
Fn = simplify(Mn*ddq - Eqn)

z  = [q; dq];
p = [m; I; l; b; g];

matlabFunction(Mn,'file','Mn','vars',{t z p});
matlabFunction(Fn,'file','Fn','vars',{t z p f});

Ml = Taylor_Matrix(Mn,z,z,[zeros(1,2*(N+1))]')

Fl = Taylor_Matrix(Fn,z,z,[zeros(1,2*(N+1))]')

A = [zeros(N+1) eye(N+1); (Ml^(-1))*jacobian(Fl,q) (Ml^(-1))*jacobian(Fl,dq)];

A = simplify(A)

B = [zeros(N+1,1); (Ml^(-1))*diff(Fl,f)];

B = simplify(B)

matlabFunction(A,'file','A_matrix','vars',{p});
matlabFunction(B,'file','B_matrix','vars',{p});

end


function M_Taylor = Taylor_Matrix(M,q,qsymb,qop)

M_Taylor = subs(M,qsymb,qop);

for i = 1:size(M,1)
    for j = 1:size(M,2)
        for k = 1:length(q)
            M_Taylor(i,j) = M_Taylor(i,j) + subs(diff(M(i,j),q(k)),qsymb,qop)*q(k);
        end
    end
end

end