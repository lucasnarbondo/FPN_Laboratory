clear all
clc

%% Set values of parameters
L12=0.0262;
R12=0.0103;
L23=0.0262;
R23=0.0358;
P4=0.42;
Q4=0;
P5=0.42;
Q5=0;
w0=1;
V0=1;
Dw=100;
Dv=20;
Dw1=100;
Dv1=20;
w01=1;
V01=1;

V_base=381;
S_base=10000;
w_base=2*pi*50;

%% Initialization
%Initialize everything as 1.0 p.u
V1=1;
V2=1;
V3=1;
w=1;
epsilon=1e-6;
delta_x=1;
k=1;
max_iter=100;

%% Main loop

while ~all(delta_x<epsilon) && k<max_iter %Loop until error<epsilon and set limit for iteration
    %Calculate the Admittance Matrix - Changes every loop (w changes)
    Z12=R12+1i*L12*w;
    Y12=1/Z12;
    Z23=R23+1i*L23*w;
    Y23=1/Z23;
    Y=[Y12 -Y12 0; -Y12 Y12+Y23 -Y23; 0 -Y23 Y23];
    
    %Calculate Power from each generator
    P1=Dw1*(w01-w);
    Q1=Dv1*(V01-abs(V1));
    P2=Dw*(w0-w);
    Q2=Dv*(V0-abs(V2));
    P3=Dw*(w0-w);
    Q3=Dv*(V0-abs(V3));
    
    %Calculate Voltage at each node using Gauss-Seidel
    V1_next=1/Y(1,1)*( ((P1-P4)-1i*(Q1-Q4))/conj(V1) - Y(1,2)*V2 -Y(1,3)*V3);
    V2_next=1/Y(2,2)*( (P2-1i*Q2)/conj(V2) - Y(2,1)*V1_next - Y(2,3)*V3);
    V3_next=1/Y(3,3)*( ((P3-P5)-1i*(Q3-Q5))/conj(V3) - Y(3,1)*V1_next - Y(3,2)*V2_next);
    
    %Calculate Voltage error
    delta_V=[abs(V1_next-V1) abs(V2_next-V2) abs(V3_next-V3)];
    
    %Set voltage for next step
    V1=V1_next;
    V2=V2_next;
    V3=V3_next;
    
    %Calculate system losses
    I12=(V1-V2)*Y12;
    I23=(V2-V3)*Y23;
    P_loss=R12*abs(I12)^2+R23*abs(I23)^2;  %Losses as R*I^2, I=DeltaV*Y
    Q_loss=w*L12*abs(I12)^2+w*L23*abs(I23)^2;  %Reactive losses for checking
    
    %Calculate next step frequency
    w_next=(2*Dw*w0+Dw1*w01-(P4+P5+P_loss))/(2*Dw+Dw1);
    
    %Calculate frequency error
    delta_w=w_next-w;   
    
    %Set next step frequency
    w=w_next;
    
    %Values to evaluate for next step
    delta_x=[delta_V abs(delta_w)];
    k=k+1;
end

I_base=S_base/V_base/sqrt(3);

%Print the results
fprintf('V1 = %.2f < %.2f\n',abs(V1)*V_base,angle(V1)*180/pi);
fprintf('V2 = %.2f < %.2f\n',abs(V2)*V_base,angle(V2)*180/pi);
fprintf('V3 = %.2f < %.2f\n',abs(V3)*V_base,angle(V3)*180/pi);
fprintf('f = %.2f\n',w*w_base/2/pi);
fprintf('S1 = %.1f%+.1fj\n',P1*S_base,Q1*S_base);
fprintf('S2 = %.1f%+.1fj\n',P2*S_base,Q2*S_base);
fprintf('S3 = %.1f%+.1fj\n',P3*S_base,Q3*S_base);
fprintf('I12 = %.2f < %.2f\n',abs(I12)*I_base,angle(I12)*180/pi);
fprintf('I23 = %.2f < %.2f\n',abs(I23)*I_base,angle(I23)*180/pi);
