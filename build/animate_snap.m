clear; clc;close all

N    = 13;
ne   = 200;
nt   = 10000;
isnap= 10;

f = fopen("OUTPUT/SEM_snapshots_V.bin","r");
u = fread(f,"float64");
u = reshape(u,nt/isnap,[]);




figure()
for i=1:nt/isnap
   plot(u(i,:));

   pause(.0009)
end

