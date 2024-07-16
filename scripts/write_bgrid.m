% Define the magnetic field on a R,Z grid at fixed toroidal angles
% Rmin, Rmax, Zmin, Zmax: Boundaries of the grid in meters
% nr, nz, nphi: Grid dimensions
% nphi here is the number of planes including 0 and 2*pi/nsym, the last (symmetry) plane is not written to file 
% The grid will be defined with resolution 
% dr = (Rmax-Rmin)/(nr-1)
% dz = (Zmax-Zmin)/(nz-1)
% dphi = 2*pi/nsym/(nphi-1)

%% Define arrays
r = linspace(Rmin,Rmax,nr);
z = linspace(Zmin,Zmax,nz);
phi = linspace(0,2*pi/nsym,nphi);

Br = zeros(nr,nz,nphi);
Bz = zeros(nr,nz,nphi);
Bphi = zeros(nr,nz,nphi);

%% Fill B arrays

%% Write files
fid = fopen([out_file,'_layout.dat'],'w');
fprintf(fid,'%d %d %d %d %f %f %f %f\n',nr,nz,nphi-1,nsym,Rmin,Rmax,Zmin,Zmax);
fclose(fid);

fid = fopen([out_file,'_r.dat'],'w');
fprintf('Writing Br\n')
write_loop2(fid,Br,nr,nz,nphi);
fclose(fid);

fid = fopen([out_file,'_z.dat'],'w');
fprintf('Writing Bz\n')
write_loop2(fid,Bz,nr,nz,nphi);
fclose(fid);

fid = fopen([out_file,'_phi.dat'],'w');
fprintf('Writing Bphi\n')
write_loop2(fid,Bphi,nr,nz,nphi);
fclose(fid);

function write_loop2(fid,B3D,nr,nz,nphi)
for i = 1:nr
    for j = 1:nz
        for k = 1:nphi-1
            fprintf(fid,'%18.12e\n',B3D(i,j,k));
        end
    end
end
end
