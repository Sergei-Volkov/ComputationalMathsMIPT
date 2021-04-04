% Волков, вариант 5
% шеститочечная параметрическая схема
tL = 0; tR = 1; tau = 1e-4;
t = tL : tau : tR;
N = length(t);

xL = 0; xR = pi/2; h = 1e-1;
x = xL : h : xR;
M = length(x);

cou = tau / h^2;
disp(cou);

[T, X] = meshgrid(t, x);
u_exact = cosh(T) .* sin(X);  % точное решение

u_explicit = solve_6dot_scheme(N, M, 0, t, tau, x, h);  % неявный метод
u_implicit = solve_6dot_scheme(N, M, 1, t, tau, x, h);  % явный метод
u_kn = solve_6dot_scheme(N, M, 0.5, t, tau, x, h);  % метод К.-Н.

figure
plot3(T, X, u_exact,'--g')
hold('on')
plot3(T, X, u_explicit,'--r')
hold('on')
plot3(T, X, u_implicit,'--b')
hold('on')
plot3(T, X, u_kn,'--y')

xlabel('t')
ylabel('x')
zlabel('u(x,t)')
grid on;


% функции
function num_solution=solve_6dot_scheme(N, M, ksi, t, tau, x, h)
    const1 = ksi*tau/h^2; const2 = (1 - ksi)*tau/h^2;
    
    num_solution = zeros(N, M);
    for m = 1:M  % задаем значения функции при t=0
        num_solution(1, m) = sin(x(m));
    end
    
    A = zeros(M);  % матрица для нахождения значения функции в следующий
    for i = 1:M    % момент времени по предыдущему
        A(i, i) = 1 + 2*const1;
        if i == 1
            A(i, i+1) = -const1;
        elseif i == M
            A(i, i-1) = -const1;
        else
            A(i, i+1) = -const1;
            A(i, i-1) = -const1;
        end
    end
    
    for n = 1:N 
        % задаем условия на производную в точках 1 и М
        num_solution(n, 2) = num_solution(n, 1) + h*cosh(t(n));
        num_solution(n, M) = num_solution(n, M-1);
        f = zeros(M, 1); % строим вектор правой части уравнения
        exp_tn = exp(t(n));
        for m = 1:M
            f(m) = f(m) + tau*exp_tn*sin(x(m));
            f(m) = f(m) + (1 - 2*const2)*num_solution(n, m);
            if m == 1
                f(m) = f(m) + const2*num_solution(n, m+1);
            elseif m == M
                f(m) = f(m) + const2*num_solution(n, m-1);
            else
                f(m) = f(m) + const2*(num_solution(n, m+1) + num_solution(n, m-1));
            end
        end
        
        if n ~= N
            num_solution(n+1, :) = A\f;  % решение СЛАУ - функция в следующий момент времени
        end
    end
end