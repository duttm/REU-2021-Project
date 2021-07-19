clear
%% parameters
% single parameter:
    % 1. k2
% load('09292019_k2_0619_t'); x = Solution(2,:); k = 0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
% load('09292019_k2_0619_t2'); x = Solution(2,:); k = 0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
% load('0214_k2.mat'); x = Solution(1,:); k=[-1.3 0 0 0 +1.3 -8]; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(1) = g; x(30:41) = fvc;% corrected detailed balance(all equal assumption); k2b = g*k2;
    % 2. k3
% load('0214_KMT.mat'); x = Solution(1,:); k=[-1 0 0 -1 +2 -6]; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(7) = f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); k3b = f*k3;
    % 3. k5
% load('0214_KMP.mat'); x = Solution(1,:); kfdq=0; x([27 28 29]) = [0 0 1];  f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(10) = 1/f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); %k5b = k5/f;
    % 5. km5
% load('0214_km5.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc(11) = f; x(30:41) = fvc; % corrected detailed balance(all equal assumption); %km5b = km5*f;
% combos:
    % 1. k1 k2
% load('0214_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
load('0211_k1k2combo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kd1 = kd2*f; k2b = g*k2;
    % 2. k2 k3
% load('0214_k2KMTcombo.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 7]) = [g f]; x(30:41) = fvc; % corrected detailed balance; k2b = g*k2; k3b = f*k3;
    % 3. k2 k5
% load('05242019_doseresp_k2k5.mat'); x = Solution(1,:); x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 10]) = [g f]; x(30:41) = fvc; k = 0; % k2*g, k5*f f>1
    % 4. k1k2k3km3k4k5km5k6
% load('12082019_all_0619t.mat'); x = Solution(4,:); fvc = ones(1,12); fvc([1 3 7:12]) = 10.^[x(1) x(26) x(27:32)]./10.^x([12 10 14:19]); 
% x(30:41) = fvc; x(42) = 0; % corrected detailed balance; k2b = g*k2;
% x([27 28 29]) = [0 0 1]; 

% % % % miscellaneous - test if needed, else archive.
% k1 k2 but nonlinear input function
% load('0228_k1k2_nonlinMg.mat'); x = Solution(1,:); k=0; x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% k1 k2, phoq mutant
% load('0321phoqt281r_initialmatrix.mat'); x = Solution(1,:); x(8) = log10(0.6); k = x(27:32);x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;
% load('0331_k1k2_steadystate'); x = Solution(1,:); x([7 8]) = [log10(25) log10(1)]; k = x(27:32);x([27 28 29]) = [0 0 1]; f = 10^x(26); g = 10^x(1); fvc = ones(1,12); fvc([1 3 6]) = [g f f]; x(30:41) = fvc; % corrected detailed balance; k1b = f*k1; kb1 = kb2/f; k2b = g*k2;

x(42) = 0; % phoPQ not inducible

% define all mutant parameters:
y = x; y(27) = 1; y(9) = -Inf; % -mgrB or mgrB constitutive
z = x; z([27 7]) = [1 0]; z(9) = -Inf; %-mgrB;-autoreg
z2 = x; z2(7) = 0; %-autoreg
z3 = x; z3(14:19) = x(14:19)+k; z4 = z3; z4(27) = 1; % z3:-phosphatase; %z4 = -phosphatase;-mgrB
% z3(7) = 0; z4(7) = 0; %deltaautoregphoqt281r

% initial vectors for +/- mgrB
X0 = zeros(1,19); X0del = X0;
%% dose response
%altered 1.7 to 2 below
mgrange = 10.^(-4.5:.05:2); % dosage range of mg
pp = zeros(length(mgrange),3); % mg = 1mM
[~, X1] =ode15s(@Edit_S_eqn, [0 15*3600], X0,{},x);
[~, X2] =ode15s(@Edit_S_eqn, [0 15*3600], X0del,{},y); 
[~, X3] =ode15s(@Edit_S_eqn, [0 15*3600], X0,{},z2);
trange = [0 15*3600];
for i = 1:length(mgrange)
    x(29) = mgrange(i); y(29) = mgrange(i); z2(29) = mgrange(i);
    [~, Xdr] = ode15s(@Edit_S_eqn,trange , X1(end,:),{},x);
    [~, Xdrdelmgr] =ode15s(@Edit_S_eqn, trange, X2(end,:),{},y);
    [~, Xdrdelautoreg] =ode15s(@Edit_S_eqn, trange, X3(end,:),{},z2);
    pp(i,:) = [Xdr(end,11)/Xdr(end,19) Xdrdelmgr(end,11)/Xdrdelmgr(end,19) Xdrdelautoreg(end,11)/Xdrdelautoreg(end,19)];
end
idx = length(mgrange);%find(mgrange>=30,1,'first');
figure(1);subplot(2,3,6); semilogx(mgrange, pp(:,1)/pp(idx,1),'r'); hold on; semilogx(mgrange, pp(:,2)/pp(idx,1),'b');
semilogx(mgrange, pp(:,3)/pp(idx,1),'k');
phopq_dr(:,1) = [0.03 0.1 0.3 1 3 10 30]; % Mg (mM)
phopq_dr(:,2) = [0.58 0.58 0.56 0.5 0.39 0.19 0.1]; % 
semilogx(phopq_dr(:,1),phopq_dr(:,2)/phopq_dr(end,2),'s')
xlabel('Signal, k_{-1} (representing [Mg^{2+}])'); ylabel('Normalized YFP:CFP');
xlim([mgrange(1) mgrange(end)])
legend('WT','\DeltamgrB','\Deltaautoreg','WT-expt')
title('Dose-response');

%% sensitivity of rrp
% % % very crude way of calculating numerical derivative
% % ddmg = (rrp(2:end,:)-rrp(1:end-1,:))./(mgrange(2:end)-mgrange(1:end-1))';
% % loggain = ddmg.*mgrange(1:end-1)'./rrp(1:end-1,:);
% % figure(3); semilogx(mgrange(1:end-1),abs(loggain)); hold on
% % xlabel('Signal, k_{-1}/k_{-1}^0'); ylabel('|\partial log[PhoP~P]/\partial log(k_{-1})|');
% % legend('Wild type','No negative feedback, MgrB = 0','No autoregulation')

%% 1. WT
% from Salazar et al: Cells were grown overnight in 2mM.
% diluted and shifted to 50mM for 2-3h
% washed and shifted to 2 or 0.01mM and imaged for 4h
% 50 mM --> 0.01 mM
x(29) = 2;
[~, X1] =ode15s(@Edit_S_eqn, [0 8*3600], X0,{},x);
x(29) = 50;
[~, X] =ode15s(@Edit_S_eqn, [0 3*3600], X1(end,:),{},x);
x(29) = 0.01;
[t, Y]=ode15s(@Edit_S_eqn, [0 20]*3600, X(end,:),{},x);
% 2mM --> 0.01 mM
x(29) = 2;
[~, X_lh2] =ode15s(@Edit_S_eqn, [0 3*3600], X1(end,:),{},x);
x(29) = 0.01;
[t_lh2, Y_lh2]=ode15s(@Edit_S_eqn, [0 20]*3600, X_lh2(end,:),{},x);
% 50mM --> 10mM
x(29) = 10;
[t_lh3, Y_lh3]=ode15s(@Edit_S_eqn, [0 20]*3600, X(end,:),{},x);
% 50mM --> 2mM
x(29) = 2;
[t_lh4, Y_lh4]=ode15s(@Edit_S_eqn, [0 20]*3600, X(end,:),{},x);
%% 2. delta mgrB, w/ autoreg
y(29) = 2;
[~, X1delmgr] =ode15s(@Edit_S_eqn, [0 8*3600], X0del,{},y);
y(29) = 50;
[~, X_delmgr] =ode15s(@Edit_S_eqn, [0 3*3600], X1delmgr(end,:),{},y);
y(29) = 0.01;
[t_delmgr, Y_delmgr]=ode15s(@Edit_S_eqn, [0 20]*3600, X_delmgr(end,:),{},y);
%% 3. delta mgrB, w/o autoreg
z(29) = 2;
[~, X1delmgrdelautoreg] =ode15s(@Edit_S_eqn, [0 8*3600], X0del,{},z);
z(29) = 50;
[~, X_delmgrdelautoreg] =ode15s(@Edit_S_eqn, [0 3*3600], X1delmgrdelautoreg(end,:),{},z);
z(29) = 0.01;
[t_delmgrdelautoreg, Y_delmgrdelautoreg]=ode15s(@Edit_S_eqn, [0 20]*3600, X_delmgrdelautoreg(end,:),{},z);
%% 4. del autoreg only
z2(29) = 2;
[~, X1delautoreg] =ode15s(@Edit_S_eqn, [0 8*3600], X0,{},z2);
z2(29) = 50;
[~, X_delautoreg] =ode15s(@Edit_S_eqn, [0 3*3600], X1delautoreg(end,:),{},z2);
z2(29) = 0.01;
[t_delautoreg, Y_delautoreg]=ode15s(@Edit_S_eqn, [0 20]*3600, X_delautoreg(end,:),{},z2);
%% 5. Constitutive expression of mgrB
y(9) = x(6)+1; %10xkbtpn2
y(29) = 2;
[~, X1const] =ode15s(@Edit_S_eqn, [0 8*3600], X0,{},y);
y(29) = 50;
[~, X_const] =ode15s(@Edit_S_eqn, [0 3*3600], X1const(end,:),{},y);
y(29) = 0.01;
[t_const, Y_const]=ode15s(@Edit_S_eqn, [0 20]*3600, X_const(end,:),{},y);
%% 7. PhoQ(T281R) - phosphatase lacking mutant
z3(17:19) = x(17:19)-5;
z4(17:19) = x(17:19)-5; z4(9) = -Inf;
z3(29) = 2;
[~, X1mono] =ode15s(@Edit_S_eqn, [0 8*3600], X0,{},z3);
z3(29) = 50;
[~, X_mono] =ode15s(@Edit_S_eqn, [0 3*3600], X1mono(end,:),{},z3);
z3(29) = 0.01;
[t_mono, Y_mono]=ode15s(@Edit_S_eqn, [0 20]*3600, X_mono(end,:),{},z3);
% PhoQ(t281R) delta mgrB
z4(29) = 2;
[~, X1monodelmgr] =ode15s(@Edit_S_eqn, [0 8*3600], X0del,{},z4);
z4(29) = 50;
[~, X_monodelmgr] =ode15s(@Edit_S_eqn, [0 3*3600], X1monodelmgr(end,:),{},z4);
z4(29) = 0.01;
[t_monodelmgr, Y_monodelmgr]=ode15s(@Edit_S_eqn, [0 20]*3600, X_monodelmgr(end,:),{},z4);

% % 1mM Mg +/-autoreg +/-PhoQ(phosphatase) barplot
% % z3(29) = 1; % phoqt281r +autoreg
% % [~, bar3] =ode15s(@phopq_0619_t, [0 12*3600], X0(end,:),{},z3);
% % z5(29) = 1; %phoqt281r -autoreg
% % [~, bar4] =ode15s(@phopq_0619_t, [0 12*3600], X0(end,:),{},z5);
% % x(29) = 1; %WT
% % [~, bar1] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},x);
% % z2(29) = 1;
% % [~, bar2] =ode15s(@phopq_0619_t, [0 12*3600], X0,{},z2);
% % 
% % figure; bar([bar1(end,11)/bar1(end,19) bar2(end,11)/bar2(end,19) bar3(end,11)/bar3(end,19) bar4(end,11)/bar4(end,19)])
%% plotting all simulations against experimental data
load('yfp data sets.mat')
yfp_3s = yfp_datasets(:,:,1);yfp_4s = yfp_datasets(:,:,2);
yfp_5s = yfp_datasets(:,:,3);yfp_6s = yfp_datasets(:,:,4);
yfp_7s = yfp_datasets(:,:,5);yfp_8s = yfp_datasets(:,:,6);
yfp_9s = yfp_datasets(:,:,7);yfp_10s = yfp_datasets(:,:,8);
yfp_11s = yfp_datasets(:,:,9);yfp_12s = yfp_datasets(:,:,10);
yfp_1s = yfp_datasets(:,:,11);yfp_2s = yfp_datasets(:,:,12);

figure(1);% PmgrB-yfp outputs
   subplot(2,3,1) % 50 --> 0.01; WT, delmgr
plot(t/60,(Y(:,11)./Y(:,19))/(Y(1,11)/Y(1,19)),'b-'); hold on;
plot(t_delmgr/60, (Y_delmgr(:,11)./Y_delmgr(:,19))/(Y(1,11)/Y(1,19)),'r-'); %Delta mgrB
legend('WT','\DeltamgrB')
    % delmgrB
plot(yfp_1s(:,1),yfp_1s(:,2)/yfp_1s(1,2),'--'); hold on
plot(yfp_3s(:,1),yfp_3s(:,2)/yfp_1s(1,2),'k--'); % expt delta mgrB
ylabel('[YFP:CFP]/[YFP:CFP]_{WT,50mM}')
xlabel('time(mins)');
xlim([0 250]); ylim([0 45])
title('50 \rightarrow 0.01 mM : P_{mgrB}')
    subplot(2,3,2) 
     % 2 --> 0.01
plot(t_lh2/60, (Y_lh2(:,11)./Y_lh2(:,19))/(Y(1,11)/Y(1,19)),'r-'); hold on;
     % 50 --> 2
plot(t_lh4/60, (Y_lh4(:,11)./Y_lh4(:,19))/(Y(1,11)/Y(1,19)),'b-'); 
     % 50 --> 10
plot(t_lh3/60, (Y_lh3(:,11)./Y_lh3(:,19))/(Y(1,11)/Y(1,19)),'k-'); 
legend('2\rightarrow0.01','50\rightarrow2','50\rightarrow10')
plot(yfp_5s(:,1),yfp_5s(:,2)/yfp_1s(1,2),'b-'); hold on
plot(yfp_6s(:,1),yfp_6s(:,2)/yfp_1s(1,2),'k-.'); hold on
plot(yfp_2s(:,1),yfp_2s(:,2)/yfp_1s(1,2),'r--'); hold on
xlabel('time(mins)');
axis([0 250 0 20])
title('WT: P_{mgrB}')
    subplot(2,3,3) % constitutive mgrB
plot(yfp_4s(:,1), yfp_4s(:,2)/yfp_1s(1,2),'--'); hold on
plot(t_const/60, (Y_const(:,11)./Y_const(:,19))/(Y(1,11)/Y(1,19)),'r-')
xlim([0 250]); ylim([0 20])
title('constitutive mgrB:P_{mgrB}')
xlabel('time(mins)');
    subplot(2,3,4) % phopq-yfp outputs
plot(t/60, (Y(:,12)./Y(:,19))/(Y(1,11)/Y(1,19)),'k-'); hold on
    % Pphopq-delmgr
plot(t_delmgr/60, (Y_delmgr(:,12)./Y_delmgr(:,19))/(Y(1,11)/Y(1,19)),'b-'); hold on
plot(yfp_9s(:,1), yfp_9s(:,2)/yfp_1s(1,2),'--')
plot(yfp_10s(:,1), yfp_10s(:,2)/yfp_1s(1,2),'r-')
ylabel('[YFP:CFP]/[YFP:CFP]_{WT,50mM}')
xlabel('time(mins)');
title('50 \rightarrow 0.01 mM : P_{phoPQ}')
xlim([0 250]); ylim([0 12])

    subplot(2,3,5)
plot(t_delmgrdelautoreg/60, (Y_delmgrdelautoreg(:,11)./Y_delmgrdelautoreg(:,19))/(Y(1,11)/Y(1,19)),'r-'); hold on % Delta mgrB 
                                                  % delta autoreg
plot(t_delautoreg/60, (Y_delautoreg(:,11)./Y_delautoreg(:,19))/(Y(1,11)/Y(1,19)),'k-'); %delta autoreg
plot(yfp_8s(:,1),yfp_8s(:,2)/yfp_1s(1,2),'r-.'); % delmgrBdelautoreg
plot(yfp_7s(:,1),yfp_7s(:,2)/yfp_1s(1,2),'k--'); % delautoreg
%legend('\DeltamgrB\Deltaautoreg','\Deltaautoreg')%,'\DeltamgrB')
ylabel('YFP:CFP (fold)'); xlabel('time (mins)');
xlim([0 250]); %ylim([0 45]);
title('\Deltaautoreg: P_{mgrB}')

% % figure(4); subplot(2,1,1) % PhoQ(T281R)
% % plot(yfp_11s(:,1), yfp_11s(:,2)/yfp_1s(1,2),'g--'); hold on
% % plot(yfp_12s(:,1), yfp_12s(:,2)/yfp_1s(1,2),'LineStyle','--','Color',[0 0.5 0]);
% % plot(t_mono/60, (Y_mono(:,11)./Y_mono(:,19))/(Y(1,11)/Y(1,19)),'g');
% % plot(t_monodelmgr/60, (Y_monodelmgr(:,11)./Y_monodelmgr(:,19))/(Y(1,11)/Y(1,19)),'Color',[0 0.5 0]);
% % xlim([0 400]);ylim([0 60]); ylabel('YFP (fold)'); xlabel('time (mins)');
% % title('PhoQ(T281R)+/-mgrB')

% % strg = 'kpdeg kmdeg ktlnA ktlnE kbtpn1 kbtpn2 f1  K1    ktlnB  k1  km1 k2  km2 k3 km3  k4 k5 km5 k6  K2     f2   kbtpn3 ktlnY kb kd k2b km2b k3b km3b k4b';
% % vars = strsplit(strg,' ');
% % figure
% % for i = 1:size(sol,2)-1
% %     if sol(end-1,i) ~= sol(end,i)
% %     subplot(5,6,i)
% %     hist(sol(1:end-2,i))
% %     set(gca,'XLim',sol([end-1 end],i))
% %     xlabel(vars{i})
% %     end
% % end

