function trim_powell_turn()
    % Initial guess [Thrust; Elevator defl; Aileron defl; Rudder defl; AoA]
    x0 = [2; 0; 0; 0; 0.1];    % [T (lbf); δe (deg); δa (deg); δr (deg); α (rad)]
    mu_fixed = 15 * pi/180;    % Bank angle µ = 15°
    tol = 1e-6;

    % Objective for turn trim
    objFun = @(x) objectiveFunctionTurn(x, mu_fixed);

    % Run Powell optimization
    [x_trim, fval] = powell_optimization(objFun, x0, tol);

    % Display results
    fprintf('Steady Turn Trim solution:\n');
    fprintf('  Thrust        = %.6f lbf\n', x_trim(1));
    fprintf('  Elevator def  = %.6f deg\n', x_trim(2));
    fprintf('  Aileron def   = %.6f deg\n', x_trim(3));
    fprintf('  Rudder def    = %.6f deg\n', x_trim(4));
    fprintf('  AoA           = %.6f rad\n', x_trim(5));
    fprintf('  Bank angle µ  = %.2f deg\n', mu_fixed*180/pi);
    fprintf('Cost = %.3e\n', fval);
end

function J = objectiveFunctionTurn(x, mu_fixed)
    % Build input vector for EOMsturn
    inVec = zeros(7,1);
    inVec(1) = x(1);      % Thrust
    inVec(2) = x(2);      % Elevator deflection
    inVec(3) = mu_fixed;  % Bank angle µ
    inVec(4) = x(3);      % Aileron deflection
    inVec(5) = x(4);      % Rudder deflection
    inVec(7) = x(5);      % Angle of attack α

    % Evaluate EOMs for turn
    x_dot = EOMsturn(inVec);

    % Cost = sum of squares of state derivatives
    J = sum(x_dot.^2);
end

function [phi_opt, fval] = golden_section_optimization(f, a, b, tol)
    gr = (sqrt(5)-1)/2;
    c = b - gr*(b - a);
    d = a + gr*(b - a);
    fc = f(c);
    fd = f(d);
    while (b - a) > tol
        if fc > fd
            a = c; c = d; fc = fd;
            d = a + gr*(b - a); fd = f(d);
        else
            b = d; d = c; fd = fc;
            c = b - gr*(b - a); fc = f(c);
        end
    end
    phi_opt = (a + b)/2;
    fval = f(phi_opt);
end

function [x, fval] = powell_optimization(f, x0, tol)
    n = numel(x0);
    D = eye(n);  x = x0;  maxIter = 50; iter = 0;
    while iter < maxIter
        x_start = x;
        for i = 1:n
            phi = @(alpha) f(x + alpha * D(:,i));
            [alpha_opt, ~] = golden_section_optimization(phi, -1, 1, tol);
            x = x + alpha_opt * D(:,i);
        end
        if norm(x - x_start) < tol
            break;
        end
        d_new = x - x_start;
        D = [D(:,2:end), d_new];
        iter = iter + 1;
    end
    fval = f(x);
end
