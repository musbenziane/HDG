clear; clc;close all

N    = 25;
ne   = 200;
nt   = 20000;
isnap= 50;

f = fopen("OUTPUT/DG_snapshots_Sigma.bin","r");
u = fread(f,"float64");
u = reshape(u,nt/isnap,[]);


figure()
for i=1:nt/isnap
    plot(u(i,:));
    pause(.005)
end