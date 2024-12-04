clc; clear all; close all;

% drone

% Exercício 1
% Implementar, testar e analisar os controladores de estabilização vertical/angular da seção 4.3 da ultima apostila do curso de Controle de Drones da Profa. Alexandra Moutinho.

% Drone parameters

%Eletric motors 
Jm = 1.4e-5*[1 1 1 1];    % Inertia
Ke = 0.0247*[1 1 1 1];    
Kt = 0.0247*[1 1 1 1];    
L  = 0.00036*[1 1 1 1];   % indutance
R  = 0.036*[1 1 1 1];     % resistence
B  = 0.00029*[1 1 1 1];
%Helices gains
KQ = 1.2191e-5*[1 1 1 1]; % Moment constant
KT = 3.5897e-5*[1 1 1 1]; % Torque constant
%Body   
J  = diag([3.8e-3,3.8e-3,7.1e-3]);            % inertia
b  = 0.24;                 % distance from motor to center of mass
m  = 2;                 % mass
%battery
Vbat = 15;

% Actuation system
% for each motor (4 motors) we send a delta (0 to 1) of Vbat to the motor
% each motor will rotate its propeller at a speed Omega_[i=1..4]_equilibrium
% each propeller will generate a force F_[i=1..4] proportional to the square of the speed
% the sum of the forces will be the total force F

% Since we need stabilization, the steady state should be at hovering, no wind
% So the only force acting on the drone is the gravity force

G_force_total = m*9.81; % total force acting on the drone
Force_each_motor = G_force_total/4; % each motor should generate this force
% Force_each_motor = KT*Omega^2
Omega_equilibrium = sqrt(Force_each_motor/KT(1)); % motor rotation speed at steady state
Omega = [1 1 1 1]'*Omega_equilibrium;
% since the propellers are placed to cancel the torque, the total torque should be zero, but for each motor+propeller:
% Q = KQ*Omega^2

% Since it is hovering, the current change should be zero and omega should be constant
% for each motor: [electrical model]
% di/dt = 1/L*(Vbat*delta - R*i - Ke*Omega) = 0
% Vbat*delta = R*i + Ke*Omega
% for each motor: [mechanical model]
% Jm*dw/dt = Kt*i - Q - B*w = 0
% i = Kt/Jm*w + Q/B
% (Q = KQ*Omega^2)
% i = Kt/Jm*w + KQ*Omega^2/B
% Vbat*delta = R*(Kt/Jm*w + KQ*Omega^2/B) + Ke*Omega
delta_equilibrium = (Ke(1)*Omega(1) + R(1)*(KQ(1)*Omega(1) + B(1))*Omega(1)/Kt(1))/Vbat;
delta = [1 1 1 1]'*delta_equilibrium;
i_equilibrium = Kt(1)/Jm(1)*Omega(1) + KQ(1)*Omega(1)^2/B(1);

% recap of linearization
%x_dot = f(x0,u0) + df/dx|{x=x0, u=u0}*(x-x0)   + df/du|{x=x0, u=u0}*(u-u0)
%x_dot - f(x0,u0) = df/dx|{x=x0, u=u0}*(x-x0)   + df/du|{x=x0, u=u0}*(u-u0)
%x_dot -|x_dot_0| =  |---------------|* x_tilde + |----------------|* u_tilde
%|---------------|=  |---------------|* x_tilde + |----------------|* u_tilde
% x_dot_tilde     =                  A* x_tilde +                  B* u_tilde

%y-y0 = dh/dx|{x=x0, u=u0}*(x-x0) + dh/du|{x=x0, u=u0}*(u-u0)
%y_tilde = C*x_tilde + D*u_tilde

% linearization of the electrical model
% I_dot = 1/L*(Vbat*delta - R*I - Ke*Omega), I_dot changes with time, as current, as delta and as Omega
I_dot_equilibrium = 1/L(1)*(Vbat*delta(1) - R(1)*i_equilibrium - Ke(1)*Omega(1));
% derivative of I_dot with respect to delta
dI_dot_ddelta = 1/L(1)*Vbat;
% derivative of I_dot with respect to Omega
dI_dot_dOmega = -Ke(1)/L(1);
% derivative of I_dot with respect to i
dI_dot_di = -R(1)/L(1);

% I_dot - I_dot_equilibrium = dI_dot_ddelta*(delta - delta_equilibrium) + dI_dot_dOmega*(Omega - Omega_equilibrium) + dI_dot_di*(i - i_equilibrium)
% delta_tilde = delta - delta_equilibrium;
% Omega_tilde = Omega - Omega_equilibrium;
% i_tilde = i - i_equilibrium;
% I_dot_tilde = I_dot - I_dot_equilibrium;
% I_dot_tilde = dI_dot_ddelta*delta_tilde + dI_dot_dOmega*Omega_tilde + dI_dot_di*i_tilde

% linearization of the mechanical model
% w_dot = 1/Jm*(Kt*i - Q - B*w), w_dot changes with time, as w, as i and as Q
w_dot_equilibrium = 1/Jm(1)*(Kt(1)*i_equilibrium - KQ(1)*Omega(1)^2 - B(1)*Omega(1));
% derivative of w_dot with respect to i
dw_dot_di = Kt(1)/Jm(1);
% derivative of w_dot with respect to Omega
dw_dot_dOmega = -2*KQ(1)*Omega(1)/Jm(1);
% derivative of w_dot with respect to Q
dw_dot_dQ = -1/Jm(1);

% w_dot - w_dot_equilibrium = dw_dot_di*(i - i_equilibrium) + dw_dot_dOmega*(Omega - Omega_equilibrium) + dw_dot_dQ*(Q - Q_equilibrium)
% i_tilde = i - i_equilibrium;
% Omega_tilde = Omega - Omega_equilibrium;
% Q_tilde = Q - Q_equilibrium;
% w_dot_tilde = w_dot - w_dot_equilibrium;
% w_dot_tilde = dw_dot_di*i_tilde + dw_dot_dOmega*Omega_tilde + dw_dot_dQ*Q_tilde

% Actuation subsystem
% (still for each motor)
% Now we can combine the electrical and mechanical linearized models [w_dot_tilde I_dot_tilde]^T = [a11 a12;a21 a22]*[w_tilde I_tilde] + [0;b2]*delta_tilde. With output [w_tilde] = [1 0]*[w_tilde I_tilde] + [0]*delta_tilde
% recap:
% w_dot_tilde = dw_dot_di*i_tilde + dw_dot_dOmega*Omega_tilde + dw_dot_dQ*Q_tilde
% I_dot_tilde = dI_dot_ddelta*delta_tilde + dI_dot_dOmega*Omega_tilde + dI_dot_di*i_tilde

% reorganizing the equations

% w_dot_tilde = dw_dot_dOmega*Omega_tilde + dw_dot_di*i_tilde + dw_dot_dQ*Q_tilde
% I_dot_tilde = dI_dot_dOmega*Omega_tilde + dI_dot_di*i_tilde + dI_dot_ddelta*delta_tilde
% for total actuation system, Q_tilde is 0, so
% [w_dot_tilde I_dot_tilde]^T = [a11 a12;a21 a22]*[w_tilde I_tilde] + [0;b2]*delta_tilde
% [w_tilde] = [1 0]*[w_tilde I_tilde] + [0]*delta_tilde
a11 = dw_dot_dOmega;
a12 = dw_dot_di;;
a21 = dI_dot_dOmega;
a22 = dI_dot_di;

b2 = dI_dot_ddelta;

% Applying the Alexandra's super hack at section "3.4.1 Subsistema de atuação", equation 30:
% (we need to prove that)
Ga_0 = (b2*a12/(a11*a22-a12*a21)); % transfer function of Actuation subsystem: G_a(s)=Omega_i(s)/Delta_i(s)
% Applying the Alexandra's super hack at section "3.4.2 Subsistema de movimento", equation 33:
% (we need to prove that)
% height, roll, pitch and yaw transfer functions are "constant/s^2"

% vertical axis subsystem
% adding every motor to make drone go up (w is vertical speed, p_D is vertical position)
%[w_dot_tilde p_D_dot_tilde]^T = [0 0;1 0]*[w_tilde p_D_tilde] + [- b_z;0]*omega_tilde_z
% y_tilde = h_tilde = [0 -1]*[w_tilde p_D_tilde]
% where omega_tilde_z subscript indicates "omega_tilde_1 + omega_tilde_2 + omega_tilde_3 + omega_tilde_4"
% but who is b_z?

%dv/dt = -w*v + R'*g_vec + F_p/m + F_a/m, but w*v on hovering, w*v = 0, R = I, F_a=[0 0 0]'' , F_p = [0 0 -sum(T_i)]', g_vec = [0 0 g]'
% T_i = KT*Omega_i^2, so -sum(T_i) = -4*KT*Omega^2
% we just need the vertical component of the acceleration, so
% dw/dt = g - 4*KT*Omega^2/m (this w is not "omega w", but the vertical speed as in "u,v,w")
% linearizing around Omega = Omega_equilibrium
% dw/dt ~ g -4*KT/m * Omega_equilibrium^2 + (-8*KT/m * Omega_equilibrium)* (Omega - Omega_equilibrium)
% dw/dt -g -4*KT/m * Omega_equilibrium^2  ~ (-8*KT/m * Omega_equilibrium)* (Omega - Omega_equilibrium)
% w_(vspeed)_dot_tilde                    = -8*KT/m * Omega_equilibrium  * Omega_tilde

% since Alexadra says
% w_(vspeed)_dot_tilde                    = -b_z * Omega_tilde, we have
b_z = (8*KT(1)/m)*Omega_equilibrium; % 0.053077
% bz is how much vertical acceleration grows if we increase the rotor speed by 1 rad/s
% so 0.053077 m/s^2 = 0.1910772 km/(h·s)
% so if we add only 1 rad/s to the rotor speed, the drone v_speed will increase
% to 0.1910772 km/h in 1 second, which seems reasonable since 1 rad/s is a really small increase
b_h = Ga_0*b_z;

Hgain_pd = [16/b_h 10/b_h]; % Proportional and derivative gains

% yaw subsystem
% extracted from Alexandra 3.2.2:
% O subsistema do eixo de guinada (yaw) corresponde às variáveis relativas a
%este eixo, nomeadamente a razão de guinada r e o ângulo de guinada ψ. A entrada
%corresponde à diferença entre as velocidades de rotação dos motores que induzem
%momento positivo da guinada (motores 2 e 4) e dos motores que induzem momento
%negativo (motores 1 e 3), Ω˜
%r = −Ω˜_1 +Ω˜_2 −Ω˜_3 +Ω˜_4
% we can have a w vector (overloading notation even more) with the drone angles:
% w_dot = [p q r]' = -J^(-1)*w X J.w + J^(-1)*Mp, where * is matrix multiplication, X is the cross product
% and Mp is the torque vector, M_p = [b*(T_4 - T_2); b*(T_1 - T_3); -Q_1 + Q_2 - Q_3 + Q_4]'
% since we are hovering,  -J^(-1)*w X J.w = 0, so w_dot = J^(-1)*Mp
% J is diagonal, and we just need the yaw component ( r_dot from w=[p_dot q_dot r_dot] ), so 
J_yaw_inv = 1/J(3,3);
%r_dot = J^(-1)* (-Q_1 + Q_2 - Q_3 + Q_4)
%Q_i=KQ*Omega_i^2, so linearizing around Omega = Omega_equilibrium for each motor
%r_dot = J_yaw_inv*(KQ(2)*Omega_equilibrium^2 - KQ(1)*Omega_equilibrium^2 - KQ(3)*Omega_equilibrium^2 + KQ(4)*Omega_equilibrium^2) +
%    J_yaw_inv*(KQ(2)*Omega_equilibrium^2*(Omega(2) - Omega_equilibrium) +
%             - KQ(1)*Omega_equilibrium^2*(Omega(1) - Omega_equilibrium) +
%             - KQ(3)*Omega_equilibrium^2*(Omega(3) - Omega_equilibrium) + 
%               KQ(4)*Omega_equilibrium^2*(Omega(4) - Omega_equilibrium)), but all KQ(i) are the same, so if we call KQ(1)*Omega_equilibrium^2 = Q_equilibrium
%r_dot - J_yaw_inv*Q_equilibrium = J_yaw_inv* Q_equilibrium*( 1 (Omega(2) - Omega_equilibrium) - 1 (Omega(1) - Omega_equilibrium) - 1 (Omega(3) - Omega_equilibrium) + 1 (Omega(4) - Omega_equilibrium))
%r_dot_tilde = J_yaw_inv* Q_equilibrium* [-1 1 -1 1]*Omega_tilde, and %r_dot_tilde = b_r *Omega_tilde:
b_r = J_yaw_inv*Omega_equilibrium; % 0.
b_psi = Ga_0*b_r;

PSIgain_pd = [ 2/b_psi 4/b_psi ]; % Proportional and derivative gains


% roll subsystem (same as pitch)
% from yaw subsystem, we have the torque vector, and we can extract the roll component
J_roll_inv = 1/J(1,1);
% p_dot = J_roll_inv*(b*(T_4 - T_2))
% T_i = KT*Omega_i^2, so linearizing around Omega = Omega_equilibrium for each motor
% p_dot = J_roll_inv*(b*(KT(4)*Omega_equilibrium^2 - KT(2)*Omega_equilibrium^2)) +
%    J_roll_inv*(b*(KT(4)*Omega_equilibrium^2*(Omega(4) - Omega_equilibrium) - KT(2)*Omega_equilibrium^2*(Omega(2) - Omega_equilibrium))), but all KT(i) are the same, so if we call KT(1)*Omega_equilibrium^2 = T_equilibrium
% p_dot - J_roll_inv*b*T_equilibrium = J_roll_inv*b*T_equilibrium*(1 (Omega(4) - Omega_equilibrium) - 1 (Omega(2) - Omega_equilibrium))
% p_dot_tilde = J_roll_inv*b*T_equilibrium*[-1 1]*Omega_tilde, and p_dot_tilde = b_p *Omega_tilde
b_p = J_roll_inv*b*Omega_equilibrium;
b_phi = Ga_0*b_p;
b_q = b_p;
b_theta = b_phi;

Delta =  [-1   -1  -1   -1;   
          0  -1   0   1;
          1   0  -1   0;
         -1   1  -1   1];


Bv= [0 0 0 0; 0 0 0 0; -b_z -b_z -b_z -b_z]
Bw = [0 -b_p 0 b_p;b_q 0 -b_q 0; -b_r b_r -b_r b_r]
B_INDI = [Bv; Bw; zeros(3,4); zeros(3,4)]
K  = inv(Delta)* pinv(B)
%u_indi = delta_equilibrium + pinv(B) *(erro)


% call simulacao2a_final
% outsimout_x will have Down
% outsimout_x1 will have reference for Down
% outsimout_x2 will will have the target for Down acceleration
% outsimout_x3 will have reference for Down acceleration
'simulation will take 20s + compile time (if needed)'
out = sim('simulacao2a_final');
'simulation is done'
close all;
figure()
% x axis is time in secods, y axis is Down, or negative height
'plot will take 10 seconds'
plot(out.tout , out.simout_x(:,1), 'LineWidth', 2); hold on;
plot(out.tout , out.simout_x1(:,1), 'LineWidth', 2); hold on;
legend('Down', 'Reference Down')
'plot is done'
title('Down (negative of height) vs Time')
xlabel('Time (s)')
ylabel('-Height (m)')
'save in png as simulacao2a_down.png'
saveas(gcf, 'simulacao2a_down.png')
'DONE'

close all;
figure()
% x axis is time in secods, y axis is acceleration
'plot will take 10 seconds'
plot(out.tout , out.simout_x2(:,1), 'LineWidth', 3); hold on;
plot(out.tout , out.simout_x3(:,1)*(-50/400), 'LineWidth', 2); hold on;
legend('Reference Down acceleration (with boost)', 'Reference Down acceleration')
'plot is done'
title('Down acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
'save in png as simulacao1a_down_acceleration.png'
saveas(gcf, 'simulacao2a_down_acceleration.png')
'DONE'


target_speed = 2;
air_constant = 0.254274793/10;
F_a = -air_constant * target_speed^2;
target_theta = asin(F_a/(m* (9.81)));


G_force_total_with_theta = m*9.81/cos(target_theta); % total force acting on the drone
Force_each_motor_with_theta = G_force_total_with_theta/4; % each motor should generate this force
Omega_equilibrium_with_theta = sqrt(Force_each_motor_with_theta/KT(1)); % motor rotation speed at steady state
Omega_with_theta = [1 1 1 1]'*Omega_equilibrium_with_theta;
delta_equilibrium_with_theta = (Ke(1)*Omega_with_theta(1) + R(1)*(KQ(1)*Omega_with_theta(1) + B(1))*Omega_with_theta(1)/Kt(1))/Vbat;

%simulacao2b
close all;
'simulation will take 20s + compile time (if needed)'
out = sim('simulacao2b');
'simulation is done'
close all;
figure()
% x axis is time in secods, y axis is Down, or negative height
'plot will take 10 seconds'
plot(out.tout , out.simout_x(:,1), 'LineWidth', 2); hold on;
plot(out.tout , out.simout_x1(:,1), 'LineWidth', 2); hold on;
legend('Down', 'Reference Down')
'plot is done'
title('Down (negative of height) vs Time')
xlabel('Time (s)')
ylabel('-Height (m)')
'save in png as simulacao2b_down.png'
saveas(gcf, 'simulacao2b_down.png')

% down acceleration is irrelevant, since the drone is now moving north
% theta accel is at simout_x3, reference theta is at simout_x5, theta is at simout_x4

close all;
figure()
% x axis is time in secods, y axis is acceleration
'plot will take 10 seconds'
plot(out.tout , out.simout_x3(:,1), 'LineWidth', 2); hold on;
legend('Reference Theta acceleration')
'plot is done'
title('Theta acceleration vs Time')
xlabel('Time (s)')
ylabel('Acceleration (rad/s^2)')
%vertical axis from -1 to 3
%axis([0 150 -1 3])
'save in png as simulacao2b_theta_acceleration.png'
saveas(gcf, 'simulacao2b_theta_acceleration.png')

close all;
figure()
% x axis is time in secods, y axis is theta
'plot will take 10 seconds'
plot(out.tout , out.simout_x4(:,1), 'LineWidth', 3); hold on;
plot(out.tout , out.simout_x5(:,1), 'LineWidth', 2); hold on;
legend('Theta', 'Reference Theta')
'plot is done'
title('Theta vs Time')
xlabel('Time (s)')
ylabel('Theta (rad)')
'save in png as simulacao2b_theta.png'
saveas(gcf, 'simulacao2b_theta.png')
'DONE'