% метод прямых
solve_lines_method();

function du=eq5(t, u)  
    xL = 0; xR = pi/2; h = 1e-1;
    x = xL : h : xR;
    M = length(x);
    
    ch_t = cosh(t);
    du = zeros(M, 1);
    for m = 1:M
        if m == 1
            du(m) = (u(m+1) - u(m))/h^2 - ch_t/h;
        elseif m == M
            du(m) = (u(m-1) - u(m))/h^2;
        else
            du(m) = (u(m+1) - 2*u(m) + u(m-1))/h^2;
        end
    end
end


function solve_lines_method()
    xL = 0; xR = pi/2; h = 1e-1;
    x = xL : h : xR;
    %disp(size(x))
    
    [t, u_lines] = ode45(@eq5, [0, 1], sin(x));
    [T, X] = meshgrid(t, x);
    
    u_exact = cosh(T) .* sin(X);
    disp(size(u_lines')); disp(size(u_exact));
    plot3(T, X, u_lines','--r')
    hold on;
    plot3(T, X, u_exact,'--g')
    
    xlabel('T')
    ylabel('X')
    zlabel('u(x,t)')
    grid on;
    
    dif = zeros(length(t));
    for i = 1:length(t)
        dif(i) = norm(u_exact(:, i) - u_lines(i, :));
    end
    
    %figure
    %plot(linspace(1,length(t_grid),length(t_grid)),dif)
    %xlabel('time layer')
    %ylabel('norm(dif)')
    %grid on;
end