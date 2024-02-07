% Dec 30 2022

% Solve Linearized 1D+1V Vlasov-Poisson system around background distribution f0(v) 
% Looking for damping/growth rate gamma=Im[w] and e-functions df(v)

% Assume f(t,x,v)=df(t,v)*exp(i*k*x) - small disturbance with just one k-vector
% Solve df_t + v*df_x - E*f0_v= Dvvvv*(f-<f>)_vvvv on [0,L]x[-Vmax,Vmax] with periodic boundary conditions in x and v,
% -k^2*phi=phi_xx=-rho=-Integral(df,-inf,inf)dv ~ exp(i*k*x) in k-space and 
% E=d(phi)/dx=i/k*rho
% So we solve
% df_t + i*k*v*df - i/k*rho*f0_v= Dvvvv*(f-<f>)_vvvv on [0,L]x[-Vmax,Vmax] with periodic boundary conditions in x (one mode in k-space only) and v

% Using Split Step 2nd Order and Exact solution for both of terms

% for Maxwell background (sigma=1) exact solutions for different k are (gamma=Im[w]=Landau damping): 
% k=0     w_true=1
% k=0.01  w_true=1.000150018759195 - 4.710297492898534*10^-2167*1i
% k=0.05  w_true=1.003761865294854 - 1.536295636089278*10^-84*1i
% k=0.1   w_true=1.015197525544101 - 2.612077823628287*10^-20*1i
% k=0.15  w_true=1.034829858605426 - 8.553101307336596*10^-9*1i
% k=0.2   w_true=1.063984340687703 - 0.00005510741476862841*1i
% k=0.25  w_true=1.105730823013256 - 0.002164078063942586*1i
% k=0.3   w_true=1.159846480591914 - 0.01262036842111716*1i
% k=0.34  w_true=1.20840237493423  - 0.02911192827881612*1i
% k=0.35  w_true=1.220953506161678 - 0.034318050858299615*1i

% Output of this function is w = Re[w]+1i*Im[w]. Note: the linearized code doesn't show a shift in Re[w] by mu*k if mu=!0

function w = Vlasov_1D_linearized_Steve_v4_Kappa(k, sigma, mu, kappa)

%close all                  
Vmax = 200;%8  % choose Vmax so that f0(Vmax) < 1e-16 (gm 10/13/23 use 50 for sigma=1, and 150 for sigma=5)      
L=2*pi/abs(k); % size of the system in x-direction
N=1;           % 2N is a number of grid points in x-direction, Linearized code has N=1
M=2^12;   % 2M is a number of grid points in v-direction 
dv=Vmax/M;
v=(-M:M-1)*dv;

dt=0.1;                           % Time step. Error in gamma is ~ dt^2
T_recurrence=L/dv; %min(L/dv,1000);
tfinal=T_recurrence*0.3;          % Stops the simulation before T_recurrence  
T_period=2*pi;%sqrt(1+4*k^2);     % We need approximate the value of one period, T _period=2pi/Re[w], for Maxwellian f0(v)=exp(-v^2/2)/sqrt(2pi)   T _period=~sqrt(1+4*k^2)
n_half_period=floor(T_period/2/dt); % Approximate half-period in terms of time steps dt 
t_first_measure=T_recurrence*0.1; % Time of "the first measurement". Make sure the transient dynamics in df(v) is gone after t_first_measure
n_first_measure=ceil(t_first_measure/dt)+1;
nsteps = round(tfinal / dt);      % number of time steps
nplot = 10000;                    % plot solution every nplot time steps

% Background distribution f0(v)
% Maxwell distribution
% mu=0;
% sigma=1;
% f0=f_gauss((v-mu)/sigma)/sigma;
% f0_v=f_gauss((v-mu)/sigma).*(-(v-mu)/sigma^3);

% Bi-Maxwellian Distribution
% weight1=0.8;          weight2=1-weight1;
% %mu1=-0.2*0.15;       mu2=4*(1+0.8*0.15);
% mu1 = 0;                mu2 = 4;
% sigma1=0.5;     sigma2=0.5;
% f0  =weight1*f_gauss((v-mu1)/sigma1)/sigma1 + weight2*f_gauss((v-mu2)/sigma1)/sigma1;
% f0_v=weight1*f_gauss((v-mu1)/sigma1).*(-(v-mu1)/sigma1^2)/sigma1 + weight2*f_gauss((v-mu2)/sigma2).*(-(v-mu2)/sigma2^2)/sigma2;

% Two-stream Distribution, v^2*gauss(v)
% mu=0;
% sigma=1;
% f0=(v-mu).^2.*f_gauss((v-mu)/sigma)/sigma;
% f0_v=(2*(v-mu) + (v-mu).^2.*(-(v-mu)/sigma^2)).*f_gauss((v-mu)/sigma)/sigma; 

% Lorentzian Distribution (slow decay f0(v) ~ 1/v^2)  Note: w=1-1i*k*sigma
% mu=0;
% sigma=1;
% f0=1/pi*sigma./((v-mu).^2+sigma^2); 
% f0_v=-1/pi*2*sigma*(v-mu)./((v-mu).^2+sigma^2).^2;

% % Kappa Distribution (slow decay f0(v) ~ 1/v^(2*(kappa+1))
% Changed on 12/13/22
% C1 = (pi*sigma^2*(kappa-1.5))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
% f0 = C1*(1+(v-mu).^2/((kappa-1.5)*sigma^2)).^(-kappa);
% f0_v = C1*2*(v-mu)*(-kappa)/((kappa-1.5)*sigma^2).*(1+(v-mu).^2/((kappa-1.5)*sigma^2)).^(-kappa-1);
% for comparison with non-integer kappa formula only!!
C1 = (pi*sigma^2*(kappa))^(-1/2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
f0 = C1*(1+(v-mu).^2/((kappa)*sigma^2)).^(-kappa-1);
f0_v = C1*2*(v-mu)*(-kappa)/((kappa)*sigma^2).*(1+(v-mu).^2/((kappa)*sigma^2)).^(-kappa-2);



% Kappa with bump
%C1=(pi*theta1^2*(kappa-0.5))^(-1/2)*exp(gammaln(kappa+1)-gammaln(kappa+0.5));
%C2=(pi*theta2^2*(kappa-0.5))^(-1/2)*exp(gammaln(kappa+1)-gammaln(kappa+0.5));
%f0=beta*C1*(1+(v-mu1).^2/(kappa-0.5)/theta1^2).^(-kappa-1) + ...
%(1-beta)*C2*(1+(v-mu2).^2/(kappa-0.5)/theta2^2).^(-kappa-1);
%f0_v=beta*C1*2*(v-mu1)*(-kappa-1)/(kappa-0.5)/theta1^2.*(1+(v-mu1).^2/(kappa-0.5)/theta1^2).^(-kappa-2) + ...
%(1-beta)*C2*2*(v-mu2)*(-kappa-1)/(kappa-0.5)/theta2^2.*(1+(v-mu2).^2/(kappa-0.5)/theta2^2).^(-kappa-2);

% Initial Condition for pertubation df(v)
% Gauss
% sigma_IC=1;%Vmax/10;
% df=f_gauss(v/sigma_IC)/sigma_IC; % gives best accuracy compared to the other options
omega_guess=1-1i/k;
df=f0_v./(omega_guess-v*k);

% Random Complex Gauss Field 
%R1=rand(1,2*M)-0.5+1i*(rand(1,2*M)-0.5); % Random vectors R1
%df=ifft(R1.*ifftshift(exp(-((-M:M-1)/(16)).^2))).*exp(-((-M:M-1)/(2/3*M)).^10); % Only ~16 harmonics in v-space  + filtering at -Vmax, Vmax

% Analytic E-function
%w_true=1.20840237493423  -1i*0.02911192827881612; % k=0.34  
%w_true=1.2209535061616785-1i*0.034318050858299615; % k=0.35 
%w_true=1-1i*8e-5;
%f_efunction=f_gauss(v).*(-v).*(1./(k.*v+w_true)+0./(k.*v-conj(w_true)));
%MAX_efunction=max(abs(f_efunction));
%df=f_efunction;

%Plot Background distribution f0(v) and I.C. df(v,t=0)
% figure(1);
% plot(v,f0);
% title(sprintf('Background distribution f0(v)'))
% 
% figure(7);
% plot(v,[real(df); imag(df); abs(df)]);
% title(sprintf('df(v,t=%g)  after %4i time steps with %ix%i grid points', 0,0,N,2*M))
% legend('Re[df]','Im[df]','Abs[df]')
% pause(1)

% Fourier Matrix for Advection Solve and Velocity cut off
FM_half_lin=exp(-1i*k*v*dt/2);
CutV=false; % Additional filtering at tails in v-space. Helps when initial df(v) does not decay to 0 at v->Vmax
CutVM=exp(-(v/(0.95*Vmax)).^32*dt/2); 
if CutV FM_half_lin=CutVM.*FM_half_lin; end
% Fourier Matrix for Hyper-viscosity term 
gamma_true=-0.03432;                   % use approximate value of gamma_true to compute Dvvvv
Dvvvv=abs(gamma_true)*(Vmax*2/M/pi)^4; % This is roughly what we need to kill recurrence for gamma=0.029
Dvvvv=0;
DvvvvFM=exp(-(ifftshift(-M:M-1)*2*pi/(2*Vmax)).^4*Dvvvv*dt); 

% Array Preallocation
t = 0;
id_max_start=0;
E_A=zeros(1,nsteps);
EE=zeros(1,nsteps);


% main time-stepping loop:
for n = 1:nsteps
    if n>nsteps break; end
    t = t + dt; 

% Advection, half time step   (Exact solution)
    df=FM_half_lin.*df;
% E-term, full time step (Exact Solution)
    if Dvvvv==0 
        rho=sum(df)*dv;
        E=1i*rho/k;
        df=df+E*f0_v*dt; 
    else
        rho=sum(df)*dv;
        E=1i*rho/k;
        df=df+E*f0_v*dt/2;    
        df=ifft(DvvvvFM.*fft(df)); 
        df=df+E*f0_v*dt/2; 
    end
% Advection, half time step   (Exact solution)    
    df=FM_half_lin.*df;
    
    EE(n)=E;
    E_A(n)=abs(E); % E(t) is measured at at t=(n-1/2)dt

    if E_A(n)/E_A(1)<1e-13 % stop the code if E_A reduced more than by 1e-13  
        nsteps=n;
        if n_first_measure>n/2 n_first_measure=round(n/2); end
    end
    
    % plot results at desired times:
    if mod(n,nplot)==0 || n==nsteps
%         figure(7);        
%         semilogy(v,abs(df),'r');
%         title(sprintf('df(v,t=%g)  after %4i time steps with %ix%i grid points', t,n,N,2*M))
% 
%         figure(8);
%         df_p=fft(df);
%         semilogy(-M:M-1,abs(fftshift(df_p)));  
%         title(sprintf('df_p(t=%g)  after %4i time steps with %ix%i grid points, T_{recurrence}=%g', t,n,N,2*M, L/dv))
        
        % Compute w=Re[w]+1i*Im[w]
        if n>n_first_measure
            % Compute gamma=Im[w]
            % Find the initial maximum for measurement of gamma
            if id_max_start==0
                T_move=n_half_period;
                [max_start, id_max_start]=max(E_A(n_first_measure-T_move:n_first_measure));
                id_max_start=n_first_measure-T_move+id_max_start-1;  % absolute index of the maximum
            end
            % Find the final maximum for measurement of gamma
            [max_finish, id_max_finish]=max(E_A(n-T_move:n));
            id_max_finish=n-T_move+id_max_finish-1; % absolute index of the maximum
            gamma=log(max_finish/max_start)/((id_max_finish-id_max_start)*dt);
            
            %Correct gamma by finding maximums more precisely 
            id_max_start_=max(id_max_start-ceil(T_move/2),1);
            [max_start_w, d_id]=max(E_A(id_max_start_:id_max_start_+T_move).*exp(-dt*((id_max_start_:id_max_start_+T_move)-1/2)*gamma));
            id_max_start_w=id_max_start_+d_id-1; % absolute index of the start maximum
            [max_finish_w, id_max_finish_w]=max(E_A(n-T_move:n).*exp(-dt*((n-T_move:n)-1/2)*gamma));
            id_max_finish_w=n-T_move+id_max_finish_w-1; % absolute index of the finish maximum
            gamma_add=log(max_finish_w/max_start_w)/((id_max_finish_w-id_max_start_w)*dt);
            gamma=gamma + gamma_add;

            % Compute Re[w] as a derivative of PHASE of E(t) between t_start=id_max_start*dt and t_finish=id_max_finish*dt
            E_PHASE=unwrap(angle(EE(1:n)));
            w_RE=(E_PHASE(id_max_finish_w)-E_PHASE(id_max_start_w))/((id_max_finish_w-id_max_start_w)*dt);

            w=abs(w_RE)+gamma*1i;
            
            % figure(5);
            % semilogy(dt*((0:n-1)+1/2),E_A(1:n),[id_max_finish-1/2 id_max_start-1/2]*dt,[max_finish max_start],'r*'); %.*exp([id_max_finish-1 id_max_start-1]*dt*gamma_true)
            % title(sprintf('E(t=%g)  after %4i time steps with %ix%i grid points    w=%1.8g+%0.8g*i', t,n,N,2*M,real(w),imag(w)))
            % 
            % % figure(6); % compensated for gamma decay/growth - has to be close to horizontal line
            % % semilogy(dt*((0:n-1)+1/2),E_A(1:n).*exp(-dt*((0:n-1)+1/2)*gamma),[id_max_finish-1/2 id_max_start-1/2]*dt,[max_finish max_start].*exp(-[id_max_finish-1/2 id_max_start-1/2]*dt*gamma),'r*',[id_max_finish_w-1/2 id_max_start_w-1/2]*dt,[max_finish_w max_start_w].*exp(-[id_max_finish_w-1/2 id_max_start_w-1/2]*dt*gamma_add),'m*');
            % % title(sprintf('E(t=%g)*exp(-gamma*t)  after %4i time steps with %ix%i grid points    w=%1.8g+%0.8g*i', t,n,N,2*M,real(w),imag(w)))
            % 
            % % figure(55);
            % plot(dt*((0:n-1)+1/2),E_PHASE);
            % title(sprintf('PHASE(E(t=%g))  after %4i time steps with %ix%i grid points    w=%1.8g+%0.8g*i', t,n,N,2*M,real(w),imag(w)))
        end
        % pause(0.01)
    end
end      

        

%--------------------------------------------------------
function out = f_gauss(v)
out = exp(-v.^2)/sqrt(pi);
return

%--------------------------------------------------------
function out = zi(z)
out=exp(-z.^2)*sqrt(pi)*1i - 2*mfun('dawson',z);  
return