function writePoten (filename, r, f, p)
n = length(r);
fid = fopen(filename, 'w');
fprintf (fid, '#\n%d\n%.16e\n%.16e\n', n, r(1), r(n));
fprintf (fid, '%.16e\t%.16e\t%.16e\n', [r f p]');
fclose(fid);