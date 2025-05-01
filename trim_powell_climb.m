function trim_powell_climb()
    % Initial guess [Thrust (lbf); Elevator deflection (deg); Angle of attack (rad)]
    x0 = [2; 0; 0.1];             % [T; delta_e; alpha]
    gamma_fixed = 5 * pi/180;     % Climb angle gamma = 5 degrees
    tol = 1e-6;

    % Objective function for climb trim
    objFun = @(x) objectiveFunctionClimb(x, gamma_fixed);

    % Run Powell optimization
    [trim_point, fval] = powell_optimization(objFun, x0, tol);

    % Display results
    fprintf('Steady Climb Trim solution:\n');
    fprintf('  Thrust   = %.6f lbf\n',        trim_point(1));
    fprintf('  Elevator = %.6f deg\n',        trim_point(2));
    fprintf('  AoA      = %.6f rad\n',        trim_point(3));
    fprintf('  Climb Angle γ = %.2f deg\n',    gamma_fixed * 180/pi);
    fprintf('Cost = %.3e\n',                   fval);
end

function J = objectiveFunctionClimb(x, gamma_fixed)
    % Build input vector: in(1)=T, in(2)=delta_e, in(3)=gamma, in(7)=alpha
    inVec = zeros(7,1);
    inVec(1) = x(1);
    inVec(2) = x(2);
    inVec(3) = gamma_fixed;
    inVec(7) = x(3);

    % Evaluate climb-specific EOMs
    x_dot = EOMsclimb(inVec);

    % Cost = sum of squares of state derivatives
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
            a = c; c = d; fc = fd;
            d = a + gr * (b - a); fd = f(d);
        else
            b = d; d = c; fd = fc;
            c = b - gr * (b - a); fc = f(c);
        end
    end

    phi_opt = (a + b) / 2;
    fval = f(phi_opt);
end

function [trim_point, fval] = powell_optimization(f, x0, tol)
    % Powell's optimization with initial directions
    n = numel(x0);
    D = eye(n);
    x = x0;
    maxIter = 50;
    iter = 0;

    while iter < maxIter
        x_start = x;
        % Line searches along each coordinate direction
        for i = 1:n
            phi = @(alpha) f(x + alpha * D(:, i));
            [alpha_opt, ~] = golden_section_optimization(phi, -1, 1, tol);
            x = x + alpha_opt * D(:, i);
        end
        % Check convergence
        if norm(x - x_start) < tol
            break;
        end
        % Update directions: drop first, add new
        d_new = x - x_start;
        D = [D(:, 2:end), d_new];
        iter = iter + 1;
    end

    trim_point = x;
    fval = f(x);
end
