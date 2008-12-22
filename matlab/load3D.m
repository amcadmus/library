function A = load3D (filename, n0, n1, n2)

fid = fopen (filename, 'r');
A = zeros(n0, n1, n2);

for i = 1:n0
    for j = 1:n1
        for k = 1:n2
            A(i,j,k) = fscanf(fid, '%f', [1]);
            %[i j k]
        end
    end
end

fclose (fid);