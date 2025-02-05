
clc; close all; clear; 

params = struct; diff_mats = struct;
%Fitting
lambda_eff = 4.5 ;%7.551;
h = 11.5 ;%12.33
% h = 11.8;
% h = 13.2;
h = 11;
h_tab =  51.5789;
rho = 1.2553e5;

one_C_Ah = 20;%Qbar = 20Ah
C_cell_nom = one_C_Ah;
C_rate = 4;
i_app_Amps = one_C_Ah*C_rate;

h_faces = 210;% 5;
% h_faces = 230;
Res_value = (3.4865-3.2393801212)/(2*i_app_Amps);  %   %in ohms - very sinsitive = Res.  LFP Internal resistance 20-50mOhms  = 0.02 to 0.05 = 2e-2 to 5e-2 to
n_cells = 42; %number of cell layers in one pouch cell
N = 24;%30;         % Number of grid in x,y-direction -  has to be even
% N = 50;

Temp_init = 23.8431;
%% Cell parameters
Ly = 150e-3; %Width whole battery
Lz = 200e-3; %length whole battery
Lx_onecell = (40+70)*1e-6;%x thickness 1 layer without ccs without sep - just anode and cathode

SOC0 = 0.3;
%v_min = 2; %2.75;
% v_max = 3.9;
sigma_cn_al_CC = 48650000 ;%3.77e7; % values from article   %chages shape of hot spot/where its located  - very sensitive
sigma_cp_copper_CC = 48650000 ;% 5.96e7;  % values from article
A_cell = Ly * Lz;


%%
Nl = N-1;
N2 = Nl * Nl;

[Dcheb_y,yy] = cheb(N); Dcheb_y = (2/Ly)*Dcheb_y; DDcheb_y = Dcheb_y^2;
[Dcheb_z,zz] = cheb(N); Dcheb_z = (2/Lz)*Dcheb_z; DDcheb_z = Dcheb_z^2;

z(1:N+1,1) = Lz/2*(zz(N+1:-1:1)+1); y(1:N+1,1) = Ly/2*(yy(N+1:-1:1)+1);
[X,Y] = meshgrid(y(2:N),z(2:N));

Lcn = 25e-6; %x thickness of current collector n
Lcp = 25e-6; %x thickness of current collector p
L_gap = 12.5e-3; %is distance from edges
L_tab_width = 48e-3;

for j = 1:N+1
    if y(j) <= L_gap
        start_tab1  = j;
    end
    if y(j) <= L_gap+L_tab_width
        start_tab2  = j;
    end
    if y(j) <= Ly-(L_gap+L_tab_width)
        start_tab3  = j;
    end
    if y(j) <= Ly-L_gap
        start_tab4  = j;
    end
end
start_tab5 = round(N);
tab_vec = zeros(N-1,1);
for j = 2:N
    if y(j) <= L_gap+L_tab_width && y(j) >= L_gap
        tab_vec(j-1) = 1;
    end
    if y(j) >= Ly-(L_gap+L_tab_width) && y(j) <= Ly-L_gap
        tab_vec(j-1) = 1;
    end
end

A_tab_1 = Lcp*(y(start_tab2)-y(start_tab1));
A_tab_2 = Lcn*(y(start_tab4)-y(start_tab3));

tab_mat = diag(tab_vec);
tab_mat_inv = eye(N-1)-diag(tab_vec);

I_sign_tabs = 1;
tab_bc_1 = I_sign_tabs .* i_app_Amps/(A_tab_1*sigma_cn_al_CC);
tab_bc_2 = I_sign_tabs .*i_app_Amps/(A_tab_2*sigma_cp_copper_CC);

%% Resistance matrix
gamma =  sigma_cn_al_CC*Lcn * sigma_cp_copper_CC*Lcp / (sigma_cn_al_CC*Lcn + sigma_cp_copper_CC*Lcp ) ;%goes in front of omega and I

Res_mat_noA = Res_value*ones(N-1,N-1);
Res_mat = Res_mat_noA*(A_cell);

%%
BC_vec_y = zeros(2,N-1);
BC_vec_z = [zeros(1,start_tab1-1),tab_bc_1*ones(1,start_tab2-start_tab1),...
    zeros(1,start_tab3-start_tab2),tab_bc_2*ones(1,start_tab4-start_tab3),...
    zeros(1,start_tab5-start_tab4);zeros(1,N-1)];

DDcheb_z_in_in = DDcheb_z(2:N,2:N);
DDcheb_z_in_ex = [DDcheb_z(2:N,1),DDcheb_z(2:N,N+1)];
Dcheb_z_ex_in =  [Dcheb_z(1,2:N);Dcheb_z(N+1,2:N)];
Dcheb_z_ex_ex = [Dcheb_z(1,1),Dcheb_z(1,N+1);Dcheb_z(N+1,1),Dcheb_z(N+1,N+1)];
DDcheb_z_ex = -DDcheb_z_in_ex*(Dcheb_z_ex_ex\Dcheb_z_ex_in);
DDcheb_z_ex_bc = DDcheb_z_in_ex*(Dcheb_z_ex_ex\BC_vec_z);

DDcheb_y_in_in = DDcheb_y(2:N,2:N); DDcheb_y_in_ex = [DDcheb_y(2:N,1),DDcheb_y(2:N,N+1)];
Dcheb_y_ex_in =  [Dcheb_y(1,2:N);Dcheb_y(N+1,2:N)];  Dcheb_y_ex_ex = [Dcheb_y(1,1),Dcheb_y(1,N+1);Dcheb_y(N+1,1),Dcheb_y(N+1,N+1)];
DDcheb_y_ex = -DDcheb_y_in_ex*(Dcheb_y_ex_ex\Dcheb_y_ex_in);
DDcheb_y_ex_bc = DDcheb_y_in_ex*(Dcheb_y_ex_ex\BC_vec_y);

eye_z = blkdiag(1,-1);
eye_y = blkdiag(1,-1);
reagain1 = blkdiag(1,0);
reagain2 = blkdiag(0,1);

Dcheb_y_ex_ex_T = eye_y*Dcheb_y_ex_ex+(h/lambda_eff)*eye(2);
Dcheb_z_ex_ex_T = eye_z*Dcheb_z_ex_ex+(h/lambda_eff)*eye(2);
DDcheb_y_ex_T_edge = -tab_mat_inv*DDcheb_y_in_ex*(Dcheb_y_ex_ex_T\(eye_y*Dcheb_y_ex_in));
DDcheb_z_ex_T = -DDcheb_z_in_ex*(Dcheb_z_ex_ex_T\(eye_z*Dcheb_z_ex_in));

Dcheb_y_ex_tab_T = eye_y*Dcheb_y_ex_ex+(h_tab/lambda_eff)*eye(2);
Dcheb_y_ex_ex_T_tab = eye_y*Dcheb_y_ex_tab_T;
DDcheb_y_ex_T_tab = -tab_mat*DDcheb_y_in_ex*((reagain1*Dcheb_y_ex_ex_T_tab+reagain2*Dcheb_y_ex_ex_T)\(eye_y*Dcheb_y_ex_in));

DDcheb_y_ex_T = DDcheb_y_ex_T_edge+ DDcheb_y_ex_T_tab ;

%% Dimensionless initial conditions
SOC_bar_s_p0 =  SOC0  .* ones(1, N2);
dUondT = -0.1e-3;

%% Functions - OCP
% U_Discharge = @(y,T) dUondT .* T ...
%     + 3.382+0.0047.*y+1.627.*exp(-81.163.*y.^1.0138)+7.6445e-8.*exp(25.36.*y.^2.469) - 8.441e-8.*exp(25.262.*y.^2.478);
V_exp_init = 3.2585; %V_exp_sudden = 3.398;
%Uocv_01  = 3.2651;
Uocv_03  =  3.3852;
shift_term = Uocv_03-V_exp_init ;
% U_Discharge = @(y,T,I_sign) dUondT .* T +0.0027*T*I_sign...
%     + 3.382+0.0047.*(1-y)+1.627.*exp(-81.163.*(1-y).^1.0138)+7.6445e-8.*exp(25.36.*(1-y).^2.469) - 8.441e-8.*exp(25.262.*(1-y).^2.478)-shift_term; %THIS HAS BEEN ADAPTED. 

U_Discharge = @(y,T,I_sign) dUondT .* T +0.0027*T*I_sign...
    + 3.382+0.0047.*(1-y)+1.627.*exp(-81.163.*(1-y).^1.0138)+7.6445e-8.*exp(25.36.*(1-y).^2.469) - 8.441e-8.*exp(25.262.*(1-y).^2.478)-shift_term; %THIS HAS BEEN ADAPTED. 

Ucathode0 =  U_Discharge(SOC0,0,1);

%% Initial conditions
T_0bar = zeros(1, N2);
V_0 = Ucathode0.* ones(1, N2) + Res_value .*i_app_Amps ;
Vrel_0 = zeros(1,N2);
ICs_0_orig = [SOC_bar_s_p0 T_0bar V_0 Vrel_0 Vrel_0];

ICs_0_orig = [T_0bar ];

%% M - mass matrix
D_pat_xrd = eye(N2, N2);
Z_pat_x = zeros(N2, N2);
M= [ D_pat_xrd  Z_pat_x   Z_pat_x
    Z_pat_x    D_pat_xrd  Z_pat_x
    Z_pat_x    Z_pat_x   Z_pat_x      ];

M = eye(N2);

T_store_100 = zeros((N-1)^2,1); T_store_500 = zeros((N-1)^2,1); T_store_1000 = zeros((N-1)^2,1); T_store_2500 = zeros((N-1)^2,1);

%% Store stuff. 
        params.N2 = N2; params.Nl= Nl; params.Res_mat= Res_mat;
        params.lambda_eff= lambda_eff;
        params.Lx_onecell= Lx_onecell;
       params.n_cells=  n_cells;
        params. h_faces=  h_faces;
        params.rho=  rho;
        I = 80;
       params.I=   I;
         % I = 80;
        params.Res_mat_noA= Res_mat_noA;
       diff_mats.DDcheb_z_in_in=  DDcheb_z_in_in;
        diff_mats.DDcheb_z_ex_T = DDcheb_z_ex_T;
        diff_mats.DDcheb_y_ex_T= DDcheb_y_ex_T;
        diff_mats.DDcheb_y_in_in= DDcheb_y_in_in;

%% Solve with ODE15s
n_loop = 50; %from Charles data switches at 50sec
% opts = odeset('Mass',M);
SOC_dim_1 = [];
T_dim_11 =[];
V_dim_1  =  []; time_set_loop = []; I_dim_1 = []; Vrelax_dim_1  =  []; Vrelax_dim_2  =  []; volt_tab_store = [];
Res_mat_vector = reshape(Res_mat,Nl^2,1);
for gg = 1:n_loop
    %     percent  = 1e2*gg/n_loop
    I_sign = (-1)^gg;
    [time_set,phi_sol] = ode15s(@(t,z)rhs(t,z,params, diff_mats),[0 n_loop],ICs_0_orig(:));
    ICs_0_orig(:) = phi_sol(end,:)';

    T_loop = phi_sol(:, 1 : N2);
   
    T_dim_11 = [T_dim_11;T_loop ];

    time_set_loop = [time_set_loop;time_set+(gg-1)*n_loop];

    if gg ==1 

    elseif gg ==2
        T_store_100 = T_loop;
    elseif gg == 10
        T_store_500 = T_loop;
    elseif gg == 20
        T_store_1000 = T_loop;
    elseif gg == 50
        T_store_2500 = T_loop;
    end
end

delta_t = max(max(T_dim_11))-min(min(T_dim_11));
%% Extract solution
% U_1 = U_Discharge(SOC_dim_1,T_dim_11);
Tn_1 = size(T_dim_11,1);
% I_model_for_exp_Amps =  I_dim_1;

%%
T_dim_1 = reshape(T_dim_11,Tn_1,N-1,N-1) + Temp_init;
t=time_set_loop;
%% Plot all three V SOC and T
f_size = 12;

%max T
nt = size(T_dim_1,1);
max_T = zeros(nt,1);
min_T = zeros(nt,1);
mean_T = zeros(nt,1);
for j = 1:nt
    max_T(j) = max(T_dim_1(j,:));
    min_T(j) = min(T_dim_1(j,:));
    mean_T(j) = mean(T_dim_1(j,:));
end

level_1 = [25,25.1,25.2,25.3,25.4,25.45,25.48,25.5,25.5,25.6,25.7,25.8,25.9];
level_2 = [28, 28.1, 28.2,28.3,28.4,28.6,29,29.2,29.4,29.5,29.6,29.8,30];
level_3 = [30.2,30.4,30.6,30.8,31,31.2,31.4,31.6,31.7,31.8,32,32.2];
level_4 = [31,31.2,31.4,31.6,31.8,32,32.2,32.4,32.6,32.8,32.9,33];

%% Plot model solution - T


fig1= figure;

hold on
T_plot_fit1 = reshape(T_store_100(end,:)+Temp_init,N-1,N-1);
[C,h] = contourf(X,Y,flip(T_plot_fit1),level_1,'showtext','on',FaceAlpha=0.5);
% contourf(X,Y,flip(T_plot_fit1),'showtext','on',FaceAlpha=0.5)
g = gca;
set(g, 'fontsize',f_size+2);
xlabel('y [m]','interpreter','latex','fontsize',f_size+5);
ylabel('z [m]','interpreter','latex','fontsize',f_size+5); hold off;
zlabel('T','interpreter','latex','fontsize',f_size);
box;
axis([0 Ly 0 Lz])
clabel(C,h, 'FontSize', 15);
colormap(autumn);
% box;

hold on
fig2 = figure;
T_plot_fit2 = reshape(T_store_500(end,:)+Temp_init,N-1,N-1);
[C,h] = contourf(X,Y,flip(T_plot_fit2),level_2,'showtext','on',FaceAlpha=0.5);
% contourf(X,Y,flip(T_plot_fit2),'showtext','on',FaceAlpha=0.5)
yticks([0 0.05 0.1 0.15 0.2]);
xticks([0 0.05 0.1 0.15]);
g = gca;
set(g, 'fontsize',f_size+2);
xlabel('y [m]','interpreter','latex','fontsize',f_size+5);
ylabel('z [m]','interpreter','latex','fontsize',f_size+5); hold off;
zlabel('T','interpreter','latex','fontsize',f_size);
% title('Temperature', sprintf('Time = %.0f s', 500),'interpreter','latex','fontsize',f_size)
%title('Temperature (t = 500s)','interpreter','latex','fontsize',f_size)
clabel(C,h, 'FontSize', 15);
colormap(autumn);
axis([0 Ly 0 Lz]);
% box;

fig3 = figure;
hold on
T_plot_fit3 = reshape(T_store_1000(end,:)+Temp_init,N-1,N-1);
[C,h] = contourf(X,Y,flip(T_plot_fit3),level_3,'showtext','on',FaceAlpha=0.5);
% contourf(X,Y,flip(T_plot_fit3),'showtext','on',FaceAlpha=0.5)
g = gca;
set(g, 'fontsize',f_size+2);
xlabel('y [m]','interpreter','latex','fontsize',f_size+5);
ylabel('z [m]','interpreter','latex','fontsize',f_size+5); hold off;
zlabel('T','interpreter','latex','fontsize',f_size);
axis([0 Ly 0 Lz])
clabel(C,h, 'FontSize', 15);
% title('Temperature', sprintf('Time = %.0f s', 1000),'interpreter','latex','fontsize',f_size)
%title('Temperature (t = 1000s)','interpreter','latex','fontsize',f_size)
colormap(autumn);
% box;

fig4 = figure;
hold on
T_plot_fit4 = reshape(T_store_2500(end,:)+Temp_init,N-1,N-1);
[C,h] = contourf(X,Y,flip(T_plot_fit4),level_4,'showtext','on',FaceAlpha=0.5);
% contourf(X,Y,flip(T_plot_fit4),'showtext','on',FaceAlpha=0.5)
g = gca;
set(g, 'fontsize',f_size+2);
xlabel('y [m]','interpreter','latex','fontsize',f_size+5);
ylabel('z [m]','interpreter','latex','fontsize',f_size+5); hold off;
zlabel('T','interpreter','latex','fontsize',f_size);
box;
axis([0 Ly 0 Lz])
% title('Temperature', sprintf('Time = %.0f s', 2500),'interpreter','latex','fontsize',f_size)
colormap(autumn);
clabel(C,h, 'FontSize', 15); 
% box;



hi = 1;
%%
    function out = rhs(~,X,params, diff_mats)
        % Rtilde = (3.4663-3.44141)/80;
        % Rtilde_2 = Rtilde/(33-23.724);
        N2 = params.N2;
        Nl = params.Nl; 
        Res_mat_relax = params.Res_mat;
        lambda_eff = params.lambda_eff;
        Lx_onecell = params.Lx_onecell;
        n_cells = params.n_cells;
         h_faces  = params. h_faces;
         rho = params.rho;
         I = params.I;
         I = 80;
        Res_mat_noA = params.Res_mat_noA;
        DDcheb_z_in_in = diff_mats.DDcheb_z_in_in;
        DDcheb_z_ex_T = diff_mats.DDcheb_z_ex_T;
        DDcheb_y_ex_T = diff_mats.DDcheb_y_ex_T;
        DDcheb_y_in_in = diff_mats.DDcheb_y_in_in;

        % Turn vectors into 2D matrix
        T = reshape(X(1 : N2,1),Nl,Nl);

        %% T
        FT_z = lambda_eff.*((DDcheb_z_in_in + DDcheb_z_ex_T)*T);
        FT_y= lambda_eff.*((DDcheb_y_in_in + DDcheb_y_ex_T)*T');

         Res_mat_set = Res_mat_noA;
        IsqR_term = Res_mat_set*(I).^2./(Lx_onecell *n_cells);
        T_dt = 1*( FT_z+ FT_y' +   IsqR_term - 1*h_faces .* T ./ n_cells ) ./rho; % really (T - Temp_init) at end


        %%
        out = [reshape(T_dt,Nl^2,1)];% reshape(FI,Nl^2,1)
        out = out(:);
    end
