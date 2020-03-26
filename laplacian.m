

function cpp = laplacian(dx,dy,c)

    Nx = size(c,1)-3;
    Ny = size(c,2)-3;

    mmx = (1/dx)^2;
    mmy = (1/dy)^2;
    
    ctmp = zeros(Nx+3,Ny+3);
    R = c(3:Nx+3,2:Ny+2);
    L = c(1:Nx+1,2:Ny+2);
    U = c(2:Nx+2,1:Ny+1);
    B = c(2:Nx+2,3:Ny+3);
    D = c(2:Nx+2,2:Ny+2);
    
    ctmp(2:Nx+2,2:Ny+2) = mmx*(R+L) + mmy*(U+B) - 2*(mmx+mmy)*D;
  
cpp = ctmp;
end