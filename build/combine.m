clear; clc;close all

N    = 13;
ne   = 200;
nt   = 20000;
isnap= 50;

f = fopen("OUTPUT/SEM_snapshots_Sigma.bin","r");
usem = fread(f,"float64");
usem = reshape(usem,nt/isnap,[]);

f = fopen("OUTPUT/DG_snapshots_Sigma.bin","r");
udg = fread(f,"float64");
udg = reshape(udg,nt/isnap,[]);

c   = [usem, udg];
txt = 'DG Domain ';
txt1= 'SEM Domain ';

figure()
for i=1:nt/isnap
    plot(c(i,:));
    xline(1201);
    x_points = [1201, 1201, 2601, 2601];  
    y_points = [-.16, .16, .16, -.16];
    color = [0, 0, 1];
    ylim([-.13 .16])
    xlim([1 2601])
    hold on;
    a = fill(x_points, y_points, color);
    a.FaceAlpha = 0.1;
    text(1800,0.12,txt)
    text(200,0.12,txt1)
    pause(.02)
    hold off
end



