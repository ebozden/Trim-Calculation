function trim_powell_level()
    % Initial guess [Thrust (lbf); Elevator deflection (deg); Angle of attack (rad)]
    x0 = [2; 0; 0.1];    % Adjust these as needed
    tol = 1e-6;

    % Run Powell optimization
    [trim_point, fval] = powell_optimization(@objectiveFunction, x0, tol);

    % Display results
    fprintf('Trim solution:\n');
    fprintf('  Thrust       = %.6f lbf\n', trim_point(1));
    fprintf('  Elevator def = %.6f deg\n', trim_point(2));
    fprintf('  AoA          = %.6f rad\n', trim_point(3));
    fprintf('Cost at trim: %.6e\n', fval);
end

function J = objectiveFunction(x)
    % Build input vector for EOMs
    inVec = zeros(7,1);
    inVec(1) = x(1);    % Thrust
    inVec(2) = x(2);    % Elevator deflection
    inVec(7) = x(3);    % Angle of attack

    % Evaluate equations of motion
    x_dot = EOMslevel(inVec);

    % Cost = sum of squares of all state derivatives
    J = sum(x_dot.^2);
end

function [phi_opt, fval] = golden_section_optimization(f, a, b, tol)
    % Golden Section search on interval [a, b]
    gr = (sqrt(5) - 1) / 2;
    c = b - gr * (b - a);
    d = a + gr * (b - a);
    fc = f(c);
    fd = f(d);

    while (b - a) > tol
        if fc > fd
            a = c;
            c = d;
            fc = fd;
            d = a + gr * (b - a);
            fd = f(d);
        else
            b = d;
            d = c;
            fd = fc;
            c = b - gr * (b - a);
            fc = f(c);
        end
    end

    phi_opt = (a + b) / 2;
    fval = f(phi_opt);
end

function [trim_point, fval] = powell_optimization(f, x0, tol)
    n = numel(x0);
    D = eye(n);              % Initial search directions
    x = x0;
    maxIter = 50;
    iter = 0;

    while iter < maxIter
        x_start = x;
        % Perform line searches along each direction
        for i = 1:n
            phi = @(alpha) f(x + alpha * D(:, i));
            [alpha_opt, ~] = golden_section_optimization(phi, -1, 1, tol);
            x = x + alpha_opt * D(:, i);
        end

        % Update directions: drop oldest, add new
        d_new = x - x_start;
        if norm(d_new) < tol || norm(x - x_start) < tol
            break;
        end
        D = [D(:, 2:end), d_new];
        iter = iter + 1;
    end

    trim_point = x;
    fval = f(x);
end
