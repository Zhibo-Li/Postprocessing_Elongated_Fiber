clc; close all; clear

alpha = 1:1:44; % the flow angles
im_i = sqrt(-1); % imaginary

the_ones = ones(1, length(alpha)); %

TheDD = [tand(alpha)-1*the_ones; -tand(alpha)]; % p, q
% TheDD = [-tand(alpha); tand(alpha)-1*the_ones; -2*tand(alpha)]; % p, q, r

figure('color', 'w'); 
set(gcf, 'Position',[100 100 800 600], 'DefaultTextInterpreter', 'latex')

for ii = 1:size(TheDD, 1)

    a = -2; b = 3; c = 0; d = TheDD(ii, :); %  coefficients of cubic equation

    k = 2; % k = 0, 1, 2 corresponding to three real roots (2 is the one we need)

    delta = sqrt((9*a*b*c*the_ones - 27*a^2*d - 2*b^3*the_ones).^2 - 4*((b^2-3*a*c)*the_ones).^3); 
    % corresponding to discriminant
    % in these cases, there IS A real root within d ~ (-1, 0)

    x = -b/(3*a)*the_ones + (2^(1/3)*(b^2-3*a*c)*the_ones)./(3*a* ...
        ((-1+sqrt(3)*im_i)/2)^k * (9*a*b*c*the_ones - 27*a^2*d - 2*b^3*the_ones + delta).^(1/3) ...
        ) + (((-1+sqrt(3)*im_i)/2)^k * (9*a*b*c*the_ones - 27*a^2*d - 2*b^3*the_ones + delta).^(1/3) ...
        )./(3*2^(1/3)*a*the_ones); % the roots

    plot(alpha, x, 'LineWidth', 2); hold on

end

xlabel('$\alpha$','FontSize', 24,'Interpreter', 'latex');
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'FontSize', 24, 'TickLabelInterpreter','latex')


