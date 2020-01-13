fid = fopen('basin.dat', 'r');
x = fread(fid, 2, 'double');
y = fread(fid, 2, 'double');
n = fread(fid, 1, 'int');
basin = fread(fid, [n,n], 'int')';
fclose(fid);
figure
imagesc(x,y,basin)
