set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot,'defaultTextInterpreter', 'latex');

close all

m0 = [1 2];
cv0 = .25;
r0 = .3;

s0 = cv0 * m0;
S0 = s0' .* [1 r0; r0 1] .* s0;

clip = @(x, l, u) min(u, max(l, x));
dist = @(m, s, r) wsdist(m0, S0, m0 + m*[m0(1) -m0(2)/2], s^2*[s0(1)^2 clip(r, -1, 1)*prod(s0); 
                                                               clip(r, -1, 1)*prod(s0) s0(2)^2]);

m = 0;
s = 1;
r = .3;

m_1 = fzero(@(m) dist(m, s, r) - .1, 0.1)*[m0(1) -m0(2)/2] + m0
m_2 = fzero(@(m) dist(m, s, r) - .2, 0.1)*[m0(1) -m0(2)/2] + m0
m_5 = fzero(@(m) dist(m, s, r) - .5, 0.1)*[m0(1) -m0(2)/2] + m0

s_1 = fzero(@(s) dist(m, s, r) - .1, 2)
s_2 = fzero(@(s) dist(m, s, r) - .2, 2)
s_5 = fzero(@(s) dist(m, s, r) - .5, 3)

S_s_1 = s_1^2*[s0(1)^2 r*prod(s0); r*prod(s0) s0(2)^2]
S_s_2 = s_2^2*[s0(1)^2 r*prod(s0); r*prod(s0) s0(2)^2]
S_s_5 = s_5^2*[s0(1)^2 r*prod(s0); r*prod(s0) s0(2)^2]

r_1 = fzero(@(r) dist(m, s, r) - .1, 0.2)
r_2 = fzero(@(r) dist(m, s, r) - .2, 0.2)
smr_5 = fzero(@(s) dist((m_2(1)+m_5(1))/2/m0(1)-1, s, r_1) - .5, 2)

S_r_1 = [s0(1)^2 r_1*prod(s0); r_1*prod(s0) s0(2)^2]
S_r_2 = [s0(1)^2 r_2*prod(s0); r_2*prod(s0) s0(2)^2]
S_smr_5 = smr_5^2*[s0(1)^2 r_2*prod(s0); r_2*prod(s0) s0(2)^2]

m_all = [m_1; m_2; m_5; repmat(m0, 5, 1); (m_2+m_5)/2]'
S_all = [repmat(S0(:), 1, 3) S_s_1(:) S_s_2(:) S_s_5(:) S_r_1(:) S_r_2(:) S_smr_5(:)];


col0 = [0.8500, 0.3250, 0.0980];
center0 = '-hexagram';
col1 = [0, 0.4470, 0.7410];
center1 = '-v';
            
figure
tiledlayout(3, 3)

xlims = [-1.1675 3.1675;
          0.0110 1.9890;
         -0.2662 2.2662];

for tile = 1:9
    m1 = m_all(:, tile)';
    S1 = reshape(S_all(:, tile), 2, 2);
    
    nexttile(tile)
    h2 = Estimator.mvncontour(m0, S0, center0, 'Color', col0);
    hold on
    h = Estimator.mvncontour(m1, S1, center1, 'Color', col1);
    hold off
    if tile == 3, ylim, end
    
    xlim([-1 3])
    ylim([-1 5])
    
    switch tile
        case 1, title('$\mathsf{W_2=0.1}$'), ylabel('\textsf{mean}')
        case 2, title('$\mathsf{W_2=0.2}$')
        case 3, title('$\mathsf{W_2=0.5}$')
        case 4, ylabel('\textsf{variance}')
        case 7, ylabel('\textsf{correlation}')
        case 9, title('\textsf{(mixed factors)}')
    end
    
end

