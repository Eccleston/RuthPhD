function [MeP1, MeP2] = bme_export(g1,g2,u1,u2,gM,upreg,upfactor)

tfinal = 86400; % Final time for simulation
species = {'SF','I0','I2','I1','M','T','M_T','Me','P1','M_P1','M_P1_T','Me_P1','P2','M_P2','M_P2_T','Me_P2','P0','M_P0','M_P0_T','Me_P0'}; % The list of all species
n = length(species);

pi.g1 = 0;
pi.g2 = 0;
pi.u1 = u1;
pi.u2 = u2;
pi.q = 21035;
pi.gM = gM;
pi.upreg = upreg;
pi.upfactor = upfactor;
%pi.upreg2 = upreg2;
x0i=zeros(n,1);
[ti,xi]=ode15s(@odes, [0 96*3600],x0i,[],pi);
% Write out the parameters
p.in1 = g1;
p.in2 = g2;
p.upreg = upreg;
p.upfactor = upfactor;
p.s1 = 1000;
p.s2 = 1000;
p.sf = 1000;

% Assign initial conditions
x0 = zeros(n,1);
x0(1) = p.sf;		% SF
x0(2) = 1.0;		% I0

% Solve the ODEs
[t,x] = ode15s(@odes,[0 tfinal],x0,[],p);

% Write out a solution structure to be returned by the function
for i = 1:n
  sol.(species{i}) = x(:,i);
end
Me_P1=x(:,12);
Me_P2=x(:,16);
Me_P0 = x(:,20);
% Produce a plot
figure;
plot(t, [Me_P1, Me_P2, Me_P0])
legend('Me-P1', 'Me-P2', 'Me-P0')

return

%%%

function dxdt = odes(t,x,p)

% Write out the parameters
in1 = p.in1;
in2 = p.in2;
upreg = p.upreg;
upfactor = p.upfactor;
s1 = p.s1;
s2 = p.s2;
sf = p.sf;

% Assign states
SF = x(1);
I0 = x(2);
I2 = x(3);
I1 = x(4);
M = x(5);
T = x(6);
M_T = x(7);
Me = x(8);
P1 = x(9);
M_P1 = x(10);
M_P1_T = x(11);
Me_P1 = x(12);
P2 = x(13);
M_P2 = x(14);
M_P2_T = x(15);
Me_P2 = x(16);
P0 = x(17);
M_P0 = x(18);
M_P0_T = x(19);
Me_P0 = x(20);

% Define reaction propensities
r_0 = (150.5 + (upreg * upfactor));
r_1 = (7.9892E-05 * M);
r_2 = 1505.0;
r_3 = (0.001725968 * T);
r_4 = ((1.662768E-09 * M) * T);
r_5 = (1.184643E-06 * M_T);
r_6 = (9.329349E-05 * Me);
r_7 = I1;
r_8 = (0.13 * P1);
r_9 = ((3.177334E-11 * P1) * M);
r_10 = (1.5856E-05 * M_P1);
r_11 = ((8.302928E-08 * P1) * M_T);
r_12 = (0.33353096 * M_P1_T);
r_13 = (0.0011091974705091 * M_P1_T);
r_14 = (0.1141804 * M_P1);
r_15 = (1.5856E-05 * Me_P1);
r_16 = I2;
r_17 = (0.13 * P2);
r_18 = ((3.177334E-11 * P2) * M);
r_19 = (3.4179E-05 * M_P2);
r_20 = ((8.302928E-08 * P2) * M_T);
r_21 = (0.718955265 * M_P2_T);
r_22 = (0.0011091974705091 * M_P2_T);
r_23 = (0.1141804 * M_P2);
r_24 = (3.4179E-05 * Me_P2);
r_25 = (10000.0 * I0);
r_26 = (0.13 * P0);
r_27 = ((3.177334E-11 * P0) * M);
r_28 = (0.0001 * M_P0);
r_29 = ((8.302928E-08 * P0) * M_T);
r_30 = (2.1035 * M_P0_T);
r_31 = (0.0011091974705091 * M_P0_T);
r_32 = (0.1141804 * M_P0);
r_33 = (0.0001 * Me_P0);

% Assign derivatives
dSF =0.0;
dI0 =0.0;
dI2 =0.0;
dI1 =0.0;
dM = r_0 - r_1 - r_4 + r_5 - r_9 + r_10 - r_18 + r_19 - r_27 + r_28;
dT = r_2 - r_3 - r_4 + r_5 + r_13 + r_22 + r_31;
dM_T = r_4 - r_5 - r_11 + r_12 - r_20 + r_21 - r_29 + r_30;
dMe = -r_6 + r_15 + r_24 + r_33;
dP1 = r_7 - r_8 - r_9 + r_10 - r_11 + r_12;
dM_P1 = r_9 - r_10 + r_13 - r_14;
dM_P1_T = r_11 - r_12 - r_13;
dMe_P1 = r_14 - r_15;
dP2 = r_16 - r_17 - r_18 + r_19 - r_20 + r_21;
dM_P2 = r_18 - r_19 + r_22 - r_23;
dM_P2_T = r_20 - r_21 - r_22;
dMe_P2 = r_23 - r_24;
dP0 = r_25 - r_26 - r_27 + r_28 - r_29 + r_30;
dM_P0 = r_27 - r_28 + r_31 - r_32;
dM_P0_T = r_29 - r_30 - r_31;
dMe_P0 = r_32 - r_33;

dxdt = [dSF; dI0; dI2; dI1; dM; dT; dM_T; dMe; dP1; dM_P1; dM_P1_T; dMe_P1; dP2; dM_P2; dM_P2_T; dMe_P2; dP0; dM_P0; dM_P0_T; dMe_P0];

return
