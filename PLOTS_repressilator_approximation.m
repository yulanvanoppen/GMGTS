%% Setup ----------------------------------------------------------------
                          
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clearvars
rng('default')


name = 'repressilator';
model =  'model_repressilator_Elowitz.txt';
system1 = ODEIQM(name, model, 'FixedParameters', ["DNAT" "kf" "Kd" "m1" "p1"]);

model =  'model_repressilator_Elowitz3.txt';
system2 = ODEIQM(name, model, 'FixedParameters', ["DNAT" "m1" "p1"]);


%% Expanded

close all

correlation = eye(6) + diag([.5 0 .5 0 .5], 1) + diag([.5 0 .5 0 .5], -1);
D = .01 * system1.k0 .* system1.k0' .* correlation;


generator = Generator(system1 ...
                      , 'N', 10 ...
                      , 't', [0:2.5:200] ...
                      , 'error_std', 1e-12 ...
                      , 'D', D ...
                      , 'observed', 1:system1.K ...
                      , 'varying', 1:system1.P ...
                      );

generator.generate()
generator.plot()

t = generator.data.t;
yexp = permute(generator.data.original, [1 3 2]);
betaexp = generator.data.beta;


%% Condensed ---------------------------------------------------------------


generator = Generator(system2 ...                                            % generator setup
                      , 'N', 10 ...
                      , 't', [0:2.5:200] ...
                      , 'error_std', 1e-12 ...
                      , 'D', D ...
                      , 'observed', 1:system2.K ...
                      , 'varying', 1:system2.P ...
                      );

generator.generate('silent', betaexp);
generator.plot

ycond = permute(generator.data.original, [1 3 2]);
betacond = generator.data.beta;


%% Plot --------------------------------------------------------------------

close all

col = parula(25);
colcond = col(1:10, :);
colexp = col(11:20, :);


figure('position', [100, 100, 700, 500])
tiledlayout(2, 3)

nexttile
h = plot(t, yexp(:, :, 6), '-');
set(h, {'color'}, num2cell(colexp, 2));
title("Complex $\boldsymbol{qd_{2,1}}$")
xlabel('Time (min)')
ylabel('Concentration (\textmu M)')
ylim([0 1])

nexttile
h = plot(t, yexp(:, :, 7), '-');
set(h, {'color'}, num2cell(colexp, 2));
title("Promoter $\boldsymbol{m_1}$")
xlabel('Time (min)')
ylabel('Concentration (\textmu M)')
ylim([0 .8])
ylim_m = ylim;

nexttile
h = plot(t, yexp(:, :, 1), '-');
set(h, {'color'}, num2cell(colexp, 2));
title("Repressor $\boldsymbol{p_1}$")
xlabel('Time (min)')
ylabel('Concentration (\textmu M)')
ylim_p = ylim;


nexttile(5)
h = plot(t, ycond(:, :, 4), '-');
set(h, {'color'}, num2cell(colcond, 2));
title("Promoter $\boldsymbol{m_1}$")
xlabel('Time (min)')
ylabel('Concentration (\textmu M)')
ylim(ylim_m)

nexttile(6)
h = plot(t, ycond(:, :, 1), '-');
set(h, {'color'}, num2cell(colcond, 2));
title("Repressor $\boldsymbol{p_1}$")
xlabel('Time (min)')
ylabel('Concentration (\textmu M)')
ylim(ylim_p)


%% TikZ
filestr = '../../writing/manuscript/figures/repressilator_approximation.tex';
matlab2tikz(filestr)

fid = fopen(filestr, 'r');
f = fread(fid,'*char')';
f = strrep(f, 'only marks, ', '');
f = strrep(f, 'xlabel style={', 'xlabel style={at={(axis description cs: 0.5, 0.05)}, ');
f = strrep(f, 'ylabel style={', 'ylabel style={at={(axis description cs: 0.115, 0.5)}, ');
f = strrep(f, 'title style={', 'title style={at={(axis description cs: 0.5, .96)}, ');
fclose(fid);
fid = fopen(filestr, 'w');
fwrite(fid, f);
fclose(fid);

