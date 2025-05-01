function x_dot = EOMslevel(in)
% Forces and Moments, and Equations of Motion (6-DOF)
% Updated: Define V_dot properly and ensure x_dot dimension match

% Preallocate output vector (10 states)
x_dot = zeros(10,1);

% Inputs
T    = in(1);      % Thrust (lbf)
de   = in(2);      % Elevator deflection (deg)
% Control surface defaults
dr   = 0;          % Rudder deflection (deg)
da   = 0;          % Aileron deflection (deg)

% Fixed flight conditions
V     = 90;        % Velocity (ft/s)
gamma = 0;         % Flight path angle (rad)
alpha = in(7);     % Angle of attack (rad)
q     = 0;         % Pitch rate (rad/s)
p     = 0;         % Roll rate (rad/s)
mu    = 0;         % Bank angle (rad)
beta  = 0;         % Sideslip angle (rad)
r     = 0;         % Yaw rate (rad/s)
chi   = 0;         % Heading angle (rad)
h     = 2000;      % Altitude (ft)

% Constants
d2r = pi/180; rho = 0.0023081; g = 32.17; m = 0.487669;
Ixx = 1.9480; Iyy = 1.5523; Izz = 1.9166;
S = 10.56;
CL0   = 0.421; Cla   = 4.59; CD0w  = 0.011;
ARw   = 7.9456; e     = 0.75; Kw    = 1/(pi*ARw*e);
Cmw   = -0.005; cgw   = -0.416; c     = 1.3333;
b     = 9.16;    lambda= 0.72955;
CL0t  = 0.76;   CD0t  = 0.002; Kt    = 0.446;
it    = 2*d2r;  Te    = 0.422;
nt    = 1;      St    = S;    cgt   = 3.5;
CLavt = 0.0969; CD0vt = 0.001; Svt   = S;
Tr    = 0.434;  nvt   = nt;   cgvt  = cgt;
Cmaf  = 0.114;  CDf   = 0.005;
Cnda  = -0.0128; Clda  = 0.244;
Clp   = -1/12 * Cla * (1 + 3*lambda)/(1 + lambda);
Clb   = -0.1;    Clr   = 0.01;

% Aerodynamic coefficients
CLw = CL0 + Cla*alpha;
CDw = CD0w + Kw*CLw^2;
E   = 2*(CL0 + Cla*alpha)/(pi*ARw);    % Downwash factor (simplified)
alphat = alpha + it + Te*de - E;
CLt = CL0t * alphat;
CDt = CD0t + Kt*CLt^2;
Cmf = Cmaf * alpha;
CDvt= CD0vt;

qb = 0.5 * rho * V^2;
Lw  = qb*S*CLw;
Dw  = qb*S*CDw;
Mw  = qb*S*c*Cmw;
Lt  = nt*qb*St*CLt;
Dt  = nt*qb*St*CDt;
Df  = qb*S*CDf;
Mf  = qb*S*c*Cmaf;
Dvt = nt*qb*Svt*CDvt;
Y   = nvt*qb*Svt*CLavt *(-beta + Tr*dr + r*cgvt/V);
Mc  = Lw*cgw*cos(alpha) + Dw*cgw*sin(alpha) + Mw ...
      - Lt*cgt*cos(alpha - q*cgt/V) ...
      - (Dt + Dvt)*cgt*sin(alpha - q*cgt/V) + Mf;
Nc  = -qb*nvt*Svt*CLavt *(-beta + Tr*dr + r*cgvt/V)*cgvt ...
      - qb*S*b*Cnda*da;
Lc  = qb*S*b^2/(2*V) * (Clp*p + 2*V/b*Clb*beta + Clr*dr + Clda*da*2*V/b);

% Equations of motion derivatives
V_dot     = 1/m *(-Dw*cos(beta) + Y*sin(beta) + T*cos(beta)*cos(alpha)) ...
            - g * sin(gamma);
gamma_dot = 1/(m*V) *(-Dw*sin(beta)*sin(mu) - Y*sin(mu)*cos(beta) ...
            + Lw*cos(mu) + T*(cos(mu)*sin(alpha) + sin(mu)*sin(beta)*cos(alpha))) ...
            - g/V * cos(gamma);
alpha_dot = q - tan(beta) * (p*cos(alpha) + r*sin(alpha)) ...
            - 1/(m*V*cos(beta)) * (Lw + T*sin(alpha)) ...
            + g*cos(gamma)*cos(mu) / (V*cos(beta));
q_dot     = Mc/Iyy + (Izz*p*r - Ixx*r*p)/Iyy;
p_dot     = Lc/Ixx + (Iyy*r*q - Izz*q*r)/Ixx;
mu_dot    = 1/cos(beta) * (p*cos(alpha) + r*sin(alpha)) ...
            + 1/(m*V) * (Dw*sin(beta)*cos(mu)*tan(gamma) ...
            + Y*tan(gamma)*cos(mu)*cos(beta) ...
            + Lw*(tan(beta) + tan(gamma)*sin(mu)) ...
            + T*(sin(alpha)*tan(gamma)*sin(mu) + sin(alpha)*tan(beta) ...
            - cos(alpha)*tan(gamma)*cos(mu)*sin(beta))) ...
            - g/V * cos(gamma)*cos(mu)*tan(beta);
beta_dot  = -r*cos(alpha) + p*sin(alpha) ...
            + 1/(m*V) * (Dw*sin(beta) + Y*cos(beta) ...
            - T*sin(beta)*cos(alpha)) ...
            + g/V * cos(gamma)*sin(mu);
r_dot     = Nc/Izz + (Ixx*p*q - Iyy*p*q)/Izz;
chi_dot   = 1/(m*V*cos(gamma)) * (Dw*sin(beta)*cos(mu) ...
            + Y*cos(mu)*cos(beta) + Lw*sin(mu) ...
            + T*(sin(mu)*sin(alpha) - cos(mu)*sin(beta)*cos(alpha)));
h_dot     = V * sin(gamma);

% Pack into output vector
x_dot(1)  = V_dot;
x_dot(2)  = gamma_dot;
x_dot(3)  = alpha_dot;
x_dot(4)  = q_dot;
x_dot(5)  = p_dot;
x_dot(6)  = mu_dot;
x_dot(7)  = beta_dot;
x_dot(8)  = r_dot;
x_dot(9)  = chi_dot;
x_dot(10) = h_dot;
end
