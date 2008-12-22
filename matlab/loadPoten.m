function [x, f, p, n, lower, upper] = loadPoten (filename)

fid = fopen (filename, 'r');
n = fscanf (fid, '#\n%d', 1);
lower = fscanf (fid, '%f', 1);
upper = fscanf (fid, '%f', 1);

A = fscanf (fid, '%f', [3 n]);
A = A';
fclose (fid);

x = A(:,1);
f = A(:,2);
p = A(:,3);
