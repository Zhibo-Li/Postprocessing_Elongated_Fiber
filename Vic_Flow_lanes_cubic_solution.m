%%%%%%%%%%%%%% Gx = Gy (only change alpha) %%%%%%%%%%%%%%%%%
clc; close all; clear

alpha = 0:1:45; % the flow angles
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



%%
%%%%%%%%%%%%%%%%% change Gx/Gy and alpha %%%%%%%%%%%%%%%%%
clc; close all; clear

alpha = 0:0.1:45; % the flow angles
Gx_Gy = 0.5:0.005:2; % the ratio of Gx to Gy
the_Gx_Gy = [1 0.895 1.4]; the_alpha = [35 42.4 33.1];
marker = {'.r','.b','.k'};
for foo = 1:3
    optim_Gx_Gy(foo) = find(Gx_Gy == the_Gx_Gy(foo)); optim_alpha(foo) = find(alpha == the_alpha(foo));
end
im_i = sqrt(-1); % imaginary

[XX, YY] = meshgrid(Gx_Gy, alpha);
the_ones = ones(length(alpha), length(Gx_Gy)); 

for_p = tand(YY).*XX-1*the_ones; % for calculating p
for_q = -tand(YY).*XX; % for calculating q

TheDD = [for_p(:), for_q(:)]; % p, q
% TheDD = [-tand(alpha); tand(alpha)-1*the_ones; -2*tand(alpha)]; % p, q, r

figure('color', 'w'); 
set(gcf, 'Position',[100 100 1000 600], 'DefaultTextInterpreter', 'latex')

for ii = 1:size(TheDD, 2)

    a = -2; b = 3; c = 0; d = TheDD(:, ii); %  coefficients of cubic equation

    k = 2; % k = 0, 1, 2 corresponding to three real roots (2 is the one we need)

    delta = sqrt((9*a*b*c*the_ones(:) - 27*a^2*d - 2*b^3*the_ones(:)).^2 - 4*((b^2-3*a*c)*the_ones(:)).^3); 
    % corresponding to discriminant
    % in these cases, there IS A real root within d ~ (-1, 0)

    x = -b/(3*a)*the_ones(:) + (2^(1/3)*(b^2-3*a*c)*the_ones(:))./(3*a* ...
        ((-1+sqrt(3)*im_i)/2)^k * (9*a*b*c*the_ones(:) - 27*a^2*d - 2*b^3*the_ones(:) + delta).^(1/3) ...
        ) + (((-1+sqrt(3)*im_i)/2)^k * (9*a*b*c*the_ones(:) - 27*a^2*d - 2*b^3*the_ones(:) + delta).^(1/3) ...
        )./(3*2^(1/3)*a*the_ones(:)); % the roots

    reshape_x = real(reshape(x, size(XX, 1), size(XX, 2)));
%     plot(YY(:,6), reshape_x(:,6), 'LineWidth', 2); hold on
    if ii == 1
        surf(XX, YY, reshape_x,'FaceColor','m','FaceAlpha',0.3,'EdgeColor','m','EdgeAlpha',0.3); hold on
        for foo = 1:3
            plot3(the_Gx_Gy(foo), the_alpha(foo), reshape_x(optim_alpha(foo), ...
                optim_Gx_Gy(foo)), marker{foo}, 'MarkerSize', 30); hold on
        end
    else
        surf(XX, YY, reshape_x,'FaceColor','c','FaceAlpha',0.3,'EdgeColor','c','EdgeAlpha',0.3); hold on
        for foo = 1:3
            plot3(the_Gx_Gy(foo), the_alpha(foo), reshape_x(optim_alpha(foo), ...
                optim_Gx_Gy(foo)), marker{foo}, 'MarkerSize', 30); hold on
        end
    end

end

legend({'p','','','','q'}, 'FontName','Times New Roman')

xlim([0.5 2]); ylim([0 45]); zlim([0 1])
xlabel('$G_x/G_y$','FontSize', 24,'Interpreter', 'latex');
ylabel('$\alpha$','FontSize', 24,'Interpreter', 'latex');
set(gca,'Box', 'On','XGrid', 'On', 'YGrid', 'On', 'FontSize', 24, 'TickLabelInterpreter','latex')
view(65.5,14.5)

f=gcf;
exportgraphics(f,'D:\Dropbox\Transfer\flow_lane_3D.png','Resolution',100)