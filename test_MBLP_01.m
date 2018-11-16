clc
clear 

A1 = [0.9985 0.009372 -0.0008435 -0.08045 0;
      -0.1503 0.9004 1.748 0.006218 0;
      0.006644 -0.02278 0.8762 -0.0002468 0;
      0.0007664 -0.002952 0.2352 1 0;
      0.01978 -0.2403 0.00612 1.934 1];
B1 = [0.002361 0.2498;
      -0.3057 0.1065;
      -0.2723 0.1415;
      -0.0348 0.01804;
      0.005046 0.003108];
A2 = [0.8846 0.2159 -6.048 -5.712 -0.06709;
      -0.1337 0.833 2.974 1.14 0.01404;
      8.402*10^-5 -1.705*10^-5 0.8488 -0.03966 2.302*10^-5;
      1.093*10^-5 -2.541*10^-6 0.2308 0.9949 3.097*10^-6;
      0.01738 -0.2358 -0.1349 1.802 0.9984];
B2 = [0.08114 -0.1287 5.859 5.569 0.06709;
      -0.01222 0.002687 -1.175 -1.122 -0.01404;
      0.004557 -0.0233 0.04707 0.03968 -2.302*10^-5;
      0.000585 -0.002991 0.006032 0.005084 -3.097*10^-6;
      0.001958 -0.003299 0.138 0.1311 0.001574];
A3 = [0.8846 0.2159 -6.048 -5.712 -0.06709;
      -0.1337 0.833 2.974 1.14 0.01404;
      8.402*10^-5 -1.705*10^-5 0.8488 -0.03966 2.302*10^-5;
      1.093*10^-5 -2.541*10^-6 0.2308 0.9949 3.097*10^-6;
      0.01738 -0.2358 -0.1349 1.802 0.9984];
B3 = [0.08114 -0.1287 5.859 5.569 0.06709;
      -0.01222 0.002687 -1.175 -1.122 -0.01404;
      0.004557 -0.0233 0.04707 0.03968 -2.302*10^-5;
      0.000585 -0.002991 0.006032 0.005084 -3.097*10^-6;
      0.001958 -0.003299 0.138 0.1311 0.001574];
Hc = 26660; % ft 

Mx = ones(5,1)*10^(5);
mx = -Mx;

Mu = ones(2,1)*10^5;
mu = -Mu;

T = 0.5; % sampling time in seconds
N = 100; %no of time steps

h_dot_des = 38; 
init_c = [0.0229 0.0413 -0.0181 0.1269]'*h_dot_des;
alpha_init = init_c(3);
theta_init = init_c(4);
u = init_c(1:2);
h_init = 24000;

Hf = 27000;
r = [0 0 0 0 Hf]';
cvx_begin
%%%Mode Dynamics%%%
variable del_hld(1,N) binary
variable del_vs(1,N) binary 
variable del_cap(1,N) binary
variable x(5,N) 
variable z11(5,N)
variable z12(2,N)
variable z21(5,N)
variable z22(5,N)
variable z13(5,N)
variable z23(5,N)

%%%%%%%%%%%%%%%%%%%
%maximize ()
subject to

x(1:5,1) == [0 alpha_init 0 theta_init h_init]';
del_hld(1,1) == 0;
del_vs(1,1) == 1;
del_cap(1,1) == 0;

for k = 1:N-1
    x(:,k+1) == A1*z11(:,k) + B1*z12(:,k) + A2*z21(:,k) + B2*z22(:,k) + A3*z13(:,k) + B3*z23(:,k);
end

%1
z11 <= repmat(Mx,1,N).*repmat(del_vs,5,1);
-z11 <= -repmat(mx,1,N).*repmat(del_vs,5,1);
z11 <= x - repmat(mx,1,N).*repmat(ones(1,N)-del_vs,5,1);
-z11 <= x + repmat(Mx,1,N).*repmat(ones(1,N)-del_vs,5,1);
%2
z12 <= repmat(Mu,1,N).*repmat(del_vs,2,1);
-z12 <= -repmat(mu,1,N).*repmat(del_vs,2,1);
z12 <= repmat(u,1,N) - repmat(mu,1,N).*repmat(ones(1,N)-del_vs,2,1);
-z12 <= repmat(u,1,N) + repmat(Mu,1,N).*repmat(ones(1,N)-del_vs,2,1);
%3
z21 <= repmat(Mx,1,N).*repmat(del_cap,5,1);
-z21 <= -repmat(mx,1,N).*repmat(del_cap,5,1);
z21 <= x - repmat(mx,1,N).*repmat(ones(1,N)-del_cap,5,1);
-z21 <= x + repmat(Mx,1,N).*repmat(ones(1,N)-del_cap,5,1);
%4
z22 <= repmat(Mx,1,N).*repmat(del_cap,5,1);
-z22 <= -repmat(mx,1,N).*repmat(del_cap,5,1);
z22 <= repmat(r,1,N) - repmat(mx,1,N).*repmat(ones(1,N)-del_cap,5,1);
-z22 <= repmat(r,1,N) + repmat(Mx,1,N).*repmat(ones(1,N)-del_cap,5,1);
%5
z13 <= repmat(Mx,1,N).*repmat(del_hld,5,1);
-z13 <= -repmat(mx,1,N).*repmat(del_hld,5,1);
z13 <= x - repmat(mx,1,N).*repmat(ones(1,N)-del_hld,5,1);
-z13 <= x + repmat(Mx,1,N).*repmat(ones(1,N)-del_hld,5,1);
%6
z23 <= repmat(Mx,1,N).*repmat(del_hld,5,1);
-z23 <= -repmat(mx,1,N).*repmat(del_hld,5,1);
z23 <= repmat(r,1,N) - repmat(mx,1,N).*repmat(ones(1,N)-del_hld,5,1);
-z23 <= repmat(r,1,N) + repmat(Mx,1,N).*repmat(ones(1,N)-del_hld,5,1);

del_cap + del_hld + del_vs == ones(1,N); % only one mode active

cvx_end

