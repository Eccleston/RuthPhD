function [MeP1,MeP2,M,T,MT] = simulateMHC_ifn(g1,g2,u1,u2,upreg)

tfinal = 10*24*3600; % Final time for simulation
species = {'P1','M','M_P1','M_T','M_P1_T','T','Me_P1','Me','P2','M_P2','M_P2_T','Me_P2'}; % The list of all species
%species = {'P1','M','M_P1','M_T','M_P1_T','T','Me_P1','Me','P2','M_P2','M_P2_T','Me_P2','P3','M_P3','M_P3_T','Me_P3'}; % The list of all species

n = length(species);

% Write out the parameters
p.g1 = g1;
p.g2 = g2;
p.u1 = u1;
p.u2 = u2;
p.q = 21035;
p.upreg = upreg;
%p.upreg2 = upreg2;
% Assign initial conditions
x0 = zeros(n,1);

% Solve the ODEs
[t,x] = ode15s(@odes_upreg,[0 tfinal],x0,[],p);

% Write out a solution structure to be returned by the function
MeP1 = x(end,7);
MeP2 = x(end,12);
M = x(end,2);
T = x(end,6);
MT = x(end,4);
return

%%%

function dxdt = odes(t,x,p)

% Write out the parameters
g1 = p.g1;
g2 = p.g2;
u1 = p.u1;
u2 = p.u2;
q = p.q;

% Assign states
P1 = x(1);
M = x(2);
M_P1 = x(3);
M_T = x(4);
M_P1_T = x(5);
T = x(6);
Me_P1 = x(7);
Me = x(8);
P2 = x(9);
M_P2 = x(10);
M_P2_T = x(11);
Me_P2 = x(12);

% Define reaction propensities
r_0 = g1;
r_1 = (0.13 * P1);
r_2 = ((3.177334E-11 * P1) * M);
r_3 = (u1 * M_P1);
r_4 = ((8.302928E-08 * P1) * M_T);
r_5 = ((u1 * q) * M_P1_T);
r_6 = (0.0011091974705091 * M_P1_T);
r_7 = (0.1141804 * M_P1);
r_8 = (u1 * Me_P1);
r_9 = g2;
r_10 = (0.13 * P2);%
r_11 =  ((3.177334E-11 * P2) * M);%((3.177334E-11 * P2) * M);
r_12 = (u2 * M_P2);
r_13 = ((8.302928E-08 * P2) * M_T);
r_14 = ((u2 * q) * M_P2_T);
r_15 = (0.0011091974705091 * M_P2_T);
r_16 = (0.1141804 * M_P2);
r_17 = (u2 * Me_P2);
r_18 = 150.5;
r_19 = (7.9892E-05 * M);
r_20 = 1505.0;
r_21 = (0.001725968 * T);
r_22 = ((1.662768E-09 * M) * T);
r_23 = (1.184643E-06 * M_T);
r_24 = (9.329349E-05 * Me);

% Assign derivatives
dP1 = r_0 - r_1 - r_2 + r_3 - r_4 + r_5;
dM = -r_2 + r_3 - r_11 + r_12 + r_18 - r_19 - r_22 + r_23;
dM_P1 = r_2 - r_3 + r_6 - r_7;
dM_T = -r_4 + r_5 - r_13 + r_14 + r_22 - r_23;
dM_P1_T = r_4 - r_5 - r_6;
dT = r_6 + r_15 + r_20 - r_21 - r_22 + r_23;
dMe_P1 = r_7 - r_8;
dMe = r_8 + r_17 - r_24;
dP2 = r_9 - r_10 - r_11 + r_12 - r_13 + r_14;
dM_P2 = r_11 - r_12 + r_15 - r_16;
dM_P2_T = r_13 - r_14 - r_15;
dMe_P2 = r_16 - r_17;

dxdt = [dP1; dM; dM_P1; dM_T; dM_P1_T; dT; dMe_P1; dMe; dP2; dM_P2; dM_P2_T; dMe_P2];

return

%%
function dxdt = odes_selfpeps(t,x,p)

% Write out the parameters
g1 = p.g1;
g2 = p.g2;
u1 = p.u1;
u2 = p.u2;
q = p.q;

g3=10000;
u3=1e-3;
% Assign states
P1 = x(1);
M = x(2);
M_P1 = x(3);
M_T = x(4);
M_P1_T = x(5);
T = x(6);
Me_P1 = x(7);
Me = x(8);
P2 = x(9);
M_P2 = x(10);
M_P2_T = x(11);
Me_P2 = x(12);
P3 = x(13); %selfpep
M_P3 = x(14);
M_P3_T = x(15);
Me_P3 = x(16);
% Define reaction propensities
r_0 = g1;
r_1 = (0.13 * P1);
r_2 = ((3.177334E-11 * P1) * M);
r_3 = (u1 * M_P1);
r_4 = ((8.302928E-08 * P1) * M_T);
r_5 = ((u1 * q) * M_P1_T);
r_6 = (0.0011091974705091 * M_P1_T);
r_7 = (0.1141804 * M_P1);
r_8 = (u1 * Me_P1);
r_9 = g2;
r_10 = (0.13 * P2);
r_11 = ((3.177334E-11 * P2) * M);
r_12 = (u2 * M_P2);
r_13 = ((8.302928E-08 * P2) * M_T);
r_14 = ((u2 * q) * M_P2_T);
r_15 = (0.0011091974705091 * M_P2_T);
r_16 = (0.1141804 * M_P2);
r_17 = (u2 * Me_P2);
r_18 = 150.5;
r_19 = (7.9892E-05 * M);
r_20 = 1505.0;
r_21 = (0.001725968 * T);
r_22 = ((1.662768E-09 * M) * T);
r_23 = (1.184643E-06 * M_T);
r_24 = (9.329349E-05 * Me);
r_25 = g3; 
r_26 = (0.13 * P3);
r_27 = ((3.177334E-11 * P3) * M);
r_28 = (u3 * M_P3);
r_29 = ((8.302928E-08 * P3) * M_T);
r_30 = ((u3 * q) * M_P3_T);
r_31 = (0.0011091974705091 * M_P3_T);
r_32 = (0.1141804 * M_P3);
r_33 = (u3 * Me_P3);
% Assign derivatives
dP1 = r_0 - r_1 - r_2 + r_3 - r_4 + r_5;
dM = -r_2 + r_3 - r_11 + r_12 + r_18 - r_19 - r_22 + r_23 - r_27 + r_28;
dM_P1 = r_2 - r_3 + r_6 - r_7;
dM_T = -r_4 + r_5 - r_13 + r_14 + r_22 - r_23 - r_29 + r_30;
dM_P1_T = r_4 - r_5 - r_6;
dT = r_6 + r_15 + r_20 - r_21 - r_22 + r_23 + r_31;
dMe_P1 = r_7 - r_8;
dMe = r_8 + r_17 - r_24 + r_33;
dP2 = r_9 - r_10 - r_11 + r_12 - r_13 + r_14;
dM_P2 = r_11 - r_12 + r_15 - r_16;
dM_P2_T = r_13 - r_14 - r_15;
dMe_P2 = r_16 - r_17;
dP3 = r_25 -r_26 -r_27 + r_28 - r_29 + r_30; % self pep
dM_P3 = r_27 - r_28 + r_31 - r_32; 
dM_P3_T = r_29 - r_30 - r_31; 
dMe_P3 = r_32 - r_33; 
dxdt = [dP1; dM; dM_P1; dM_T; dM_P1_T; dT; dMe_P1; dMe; dP2; dM_P2; dM_P2_T; dMe_P2; dP3; dM_P3; dM_P3_T; dMe_P3];

return

function dxdt = odes_upreg(t,x,p)

% Write out the parameters
g1 = p.g1;
g2 = p.g2;
u1 = p.u1;
u2 = p.u2;
q = p.q;
upreg = p.upreg;
%upreg2 = p.upreg2;
% Assign states
P1 = x(1);
M = x(2);
M_P1 = x(3);
M_T = x(4);
M_P1_T = x(5);
T = x(6);
Me_P1 = x(7);
Me = x(8);
P2 = x(9);
M_P2 = x(10);
M_P2_T = x(11);
Me_P2 = x(12);

% Define reaction propensities
r_0 = g1;
r_1 = (0.13 * P1);
r_2 = ((3.177334E-11 * P1) * M);
r_3 = (u1 * M_P1);
r_4 = ((8.302928E-08 * P1) * M_T);
r_5 = ((u1 * q) * M_P1_T);
r_6 = (0.0011091974705091 * M_P1_T);
r_7 = (0.1141804 * M_P1);
r_8 = (u1 * Me_P1);
r_9 = g2;
r_10 = (0.13 * P2);%
r_11 =  ((3.177334E-11 * P2) * M);%((3.177334E-11 * P2) * M);
r_12 = (u2 * M_P2);
r_13 = ((8.302928E-08 * P2) * M_T);
r_14 = ((u2 * q) * M_P2_T);
r_15 = (0.0011091974705091 * M_P2_T);
r_16 = (0.1141804 * M_P2);
r_17 = (u2 * Me_P2);
r_18 = upreg*150.5;
r_19 = (7.9892E-05 * M);
r_20 = upreg*1505.0;
r_21 = (0.001725968 * T);
r_22 = upreg*((1.662768E-09 * M) * T);
r_23 = (1.184643E-06 * M_T);
r_24 = (9.329349E-05 * Me);

% Assign derivatives
dP1 = r_0 - r_1 - r_2 + r_3 - r_4 + r_5;
dM = -r_2 + r_3 - r_11 + r_12 + r_18 - r_19 - r_22 + r_23;
dM_P1 = r_2 - r_3 + r_6 - r_7;
dM_T = -r_4 + r_5 - r_13 + r_14 + r_22 - r_23;
dM_P1_T = r_4 - r_5 - r_6;
dT = r_6 + r_15 + r_20 - r_21 - r_22 + r_23;
dMe_P1 = r_7 - r_8;
dMe = r_8 + r_17 - r_24;
dP2 = r_9 - r_10 - r_11 + r_12 - r_13 + r_14;
dM_P2 = r_11 - r_12 + r_15 - r_16;
dM_P2_T = r_13 - r_14 - r_15;
dMe_P2 = r_16 - r_17;

dxdt = [dP1; dM; dM_P1; dM_T; dM_P1_T; dT; dMe_P1; dMe; dP2; dM_P2; dM_P2_T; dMe_P2];

return

%%