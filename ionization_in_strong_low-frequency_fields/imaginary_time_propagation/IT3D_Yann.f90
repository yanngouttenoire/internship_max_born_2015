Program IT3D

! This F90 code uses the imaginery time propagation technique to calculate the ground state of a 3D potential.
! The propagation is performed using the split operator method (time step dtfs).
! The binding potential is a harmonic potential for one electron with harmonic frequencies wx, wy and wz.
! Everything is in atomic units except when specified differently.
! Output : ground state wave function and ground state energy.
! The code is using fftw-3 (http://www.fftw.org) for the evaluation of Fourier transforms.
!
! Eric Charron, Universite Paris-Sud, CORINF, 2015.

!!
!! DECLARATIONS
!!

implicit none

integer :: nx,ny,nz,k,itmax,it,it_ch,kx,ky,kz,nsize

double precision :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,dtfs,dt,rm,eps,en,en0,variation,wx,wy,wz

double precision, parameter :: twopi=2.0d0*3.141592653589793d0
double precision, parameter :: conv_au_fs=2.418884326445171d-2,conv_au_eV=27.2113834496264d0

double precision, dimension(:), allocatable :: x,y,z,rkx,rky,rkz

double precision, dimension(:,:,:), allocatable :: V,rk3D

double complex :: ci

double complex, dimension(:), allocatable :: akx,aky,akz

double complex, dimension(:,:,:), allocatable :: ak3D,Vprop,psi,psi0

logical :: cont

integer, Parameter :: FFTW_FORWARD=-1,FFTW_BACKWARD=1
integer, Parameter :: FFTW_ESTIMATE=0,FFTW_MEASURE=1

integer*8 :: planf,planb

real, dimension(4) :: charge   
real, dimension(3) :: bondLength

!!
!! BEGIN
!!

write(*,*)
write(*,*) 'Running...'
write(*,*)

!Input parameters
nx=64           !Number of grid points for the electrons
ny=64           !Number of grid points for the electrons
nz=64           !Number of grid points for the electrons
xmax=5.0d0      !Grid maximum (au)
ymax=5.0d0      !Grid maximum (au)
zmax=5.0d0      !Grid maximum (au)
dtfs=0.01d0     !Time step (fs)
rm=1.0d0        !Electron mass (au)
wx=2.0d0        !Harmonic frequency (au)
wy=1.0d0        !Harmonic frequency (au)
wz=3.0d0        !Harmonic frequency (au)

!Conversion
dt=dtfs/conv_au_fs !Time step (au)

!Define grid
Allocate(x(nx),y(ny),z(nz))
xmin=-xmax !x-grid minimum (au)
ymin=-ymax !y-grid minimum (au)
zmin=-zmax !z-grid minimum (au)
dx=(xmax-xmin)/dble(nx-1)
dy=(ymax-ymin)/dble(ny-1)
dz=(zmax-zmin)/dble(nz-1)
do k=1,nx
 x(k)=xmin+dble(k-1)*dx
enddo
do k=1,ny
 y(k)=ymin+dble(k-1)*dy
enddo
do k=1,nz
 z(k)=zmin+dble(k-1)*dz
enddo

!Number of points
nsize=nx*ny*nz

!Define Potential

 charge(1)=0.2
 charge(2)=0.3
 charge(3)=0.3
 charge(4)=0.2
!Covalent bond length

 bondLength(1)=2.0
 bondLength(2)=1.0
 bondLength(3)=3.0
 !bondLength(1)=0.
 !bondLength(2)=0.
 !bondLength(3)=0.

Allocate(V(nx,ny,nz))
open(1,file='potential.dat',status='unknown') !Saving the potential in eV
do kx=1,nx
 do ky=1,ny
  do kz=1,nz
   V(kx,ky,kz) = -charge(1)/sqrt(x(kx)*x(kx)+y(ky)*y(ky)+z(kz)*z(kz))&
-charge(2)/sqrt((x(kx)-bondLength(1))*(x(kx)-bondLength(1))+y(ky)*y(ky)+z(kz)*z(kz))&
-charge(3)/sqrt(x(kx)*x(kx)+(y(ky)-bondLength(2))*(y(ky)-bondLength(2))+z(kz)*z(kz))&
-charge(4)/sqrt(x(kx)*x(kx)+y(ky)*y(ky)+(z(kz)-bondLength(3))*(z(kz)-bondLength(3)))
   write(1,fmt="(4(e13.6,1x))") x(kx),y(ky),z(kz),V(kx,ky,kz)*conv_au_eV
  enddo
 enddo
enddo
close(1)

!Calculate ground electronic state using imaginary time propagation

Allocate(akx(nx),aky(ny),akz(nz),ak3D(nx,ny,nz),Vprop(nx,ny,nz),psi(nx,ny,nz),rkx(nx),rky(ny),rkz(nz),rk3D(nx,ny,nz))
ci=dcmplx(1.0d0,0.0d0) !imaginary time propagation
eps=1.0d-7             !convergence criterion
itmax=10000000         !maximum number of time iterations allowed

!Initialize fft subroutines
call dfftw_plan_dft_3d(planf,nx,ny,nz,psi,psi,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_plan_dft_3d(planb,nx,ny,nz,psi,psi,FFTW_BACKWARD,FFTW_ESTIMATE)

call Defkin_1D(nx,dx,rm,rkx,akx) !Prepare the z-kinetic propagator
call Defkin_1D(ny,dy,rm,rky,aky) !Prepare the z-kinetic propagator
call Defkin_1D(nz,dz,rm,rkz,akz) !Prepare the z-kinetic propagator

call Defpot_3D(dt/2.0d0,Vprop)        !Prepare the potential propagator

do kx=1,nx
 do ky=1,ny
  do kz=1,nz
   rk3D(kx,ky,kz)=rkx(kx)+rky(ky)+rkz(kz)
   ak3D(kx,ky,kz)=akx(kx)*aky(ky)*akz(kz)/dble(nsize)      !FFT Normalization
  enddo
 enddo
enddo

!Initial guess for the wave function = constant
psi=dcmplx(1.0d0,0.0d0)

!Normalize
call norm(psi)

!Calculate initial energy
call Energy3D(psi,en)
write(*,fmt="(' iteration = ',i5,' - Energy = 'e20.13,' au.')") 0,en

!Save initial wave function
Allocate(psi0(nx,ny,nz))
psi0=psi

!loop of time propagation
it=0
it_ch=0
cont=.true.

do while((it < itmax).and.(cont))

 it=it+1

 !Apply 1/2-the potential propagator
 call proppot_3D(psi)

 !Apply the kinetic propagator
 call propkin_3D(psi)

 !Apply 1/2-the potential propagator
 call proppot_3D(psi)

 !Normalize
 call norm(psi)

 !Calculate the energy
 en0=en
 call Energy3D(psi,en)
 variation=abs((en-en0)/en)
 write(*,fmt="(' iteration = ',i5,' - Energy = 'e20.13,' au - Variation = ',1pe9.2,'.')") it,en,variation

 if (en > en0) then

  !Decrease time step for better convergence
  dt=dt/1.1d0
  write(*,*) 'Warning: energy increases : decreasing the time step for better convergence.'
  write(*,*) 'Now using dt = ',dt*conv_au_fs,' fs.'

  !Restore previous wave function
  psi=psi0
  en=en0

  !Recalculate propagators
  call Defkin_1D(nx,dx,rm,rkx,akx) !Prepare the z-kinetic propagator
  call Defkin_1D(ny,dy,rm,rky,aky) !Prepare the z-kinetic propagator
  call Defkin_1D(nz,dz,rm,rkz,akz) !Prepare the z-kinetic propagator
  call Defpot_3D(dt/2.0d0,Vprop)   !Prepare the potential propagator
  do kx=1,nx
   do ky=1,ny
    do kz=1,nz
     rk3D(kx,ky,kz)=rkx(kx)+rky(ky)+rkz(kz)
     ak3D(kx,ky,kz)=akx(kx)*aky(ky)*akz(kz)/dble(nsize)      !FFT Normalization
    enddo
   enddo
  enddo

 else

  !Save current wave function
  psi0=psi

 endif

 if (variation < eps) then
 
  !The calculation has converged with this specific time step
  !We now try to decrease the time step for better convergence

  if (it > it_ch+1) then

   !This is done only if the calculation is not fully converged
   it_ch=it
   dt=dt/1.1d0
   write(*,*) 'Converged with the specified time-step.'
   write(*,*) 'Decreasing the time step for better convergence.'
   write(*,*) 'Now using dt = ',dt*conv_au_fs,' fs.'

   !Calculate propagators
   call Defkin_1D(nx,dx,rm,rkx,akx) !Prepare the z-kinetic propagator
   call Defkin_1D(ny,dy,rm,rky,aky) !Prepare the z-kinetic propagator
   call Defkin_1D(nz,dz,rm,rkz,akz) !Prepare the z-kinetic propagator
   call Defpot_3D(dt/2.0d0,Vprop)   !Prepare the potential propagator
   do kx=1,nx
    do ky=1,ny
     do kz=1,nz
      rk3D(kx,ky,kz)=rkx(kx)+rky(ky)+rkz(kz)
      ak3D(kx,ky,kz)=akx(kx)*aky(ky)*akz(kz)/dble(nsize)      !FFT Normalization
     enddo
    enddo
   enddo

  else

   cont=.false.
   write(*,*)
   write(*,fmt="('The calculation of the ground state is now fully converged.')")
   write(*,fmt="('Found ground state in ',i5,' iterations - Energy = ',e15.8,' au.')") it,en
   write(*,*)

   !Transforming to real-valued eigenstate
   psi(:,:,:)=dcmplx(real(psi(:,:,:)),0.0d0)
   call norm(psi) !Normalization

   open(1,file='ground_state_wf.dat',status='unknown')
   do kx=1,nx
    do ky=1,ny
     do kz=1,nz
      write(1,fmt="(3(e13.6,1x),e23.16)") x(kx),y(ky),z(kz),real(psi(kx,ky,kz))
     enddo
    enddo
   enddo
   close(1)

  endif

 endif

enddo

if (cont) then
 write(*,fmt="('Failed...')")
 write(*,*)
endif

Contains

!--- SUBROUTINE ---!

   Subroutine Defkin_1D(ndim,step,mass,rks1d,aks1d)

   implicit none

   integer, intent(in) :: ndim

   double complex, dimension(:), intent(out) :: aks1d

   integer :: i

   double precision :: cons

   double precision, intent(in) :: step,mass

   double precision, dimension(:), intent(out) :: rks1d

   cons=twopi/(step*dble(ndim)) !grid step in k-space

   do i=1,ndim
    !calculate momentum
    if (i <= (ndim/2-1)) then
     rks1d(i)=dble(i-1)*cons
    else
     rks1d(i)=-dble(ndim+1-i)*cons
    endif
    !calculate (k**2)/(2m)
    rks1d(i)=rks1d(i)*rks1d(i)/(2.0d0*mass)
    !calculate 1D kinetik propagator 
    aks1d(i)=exp(-ci*rks1d(i)*dt)
   enddo

   End Subroutine Defkin_1D

!--- SUBROUTINE ---!

   Subroutine Defpot_3D(dt,Vprop)

   Implicit none

   double complex, dimension(:,:,:), intent(out) :: Vprop

   double precision, intent(in) :: dt

   !Define potential propagator
   Vprop(:,:,:)=exp(-ci*dt*V(:,:,:))

   End Subroutine Defpot_3D

!--- SUBROUTINE ---!

   Subroutine norm(psi)

   Implicit none

   double complex, dimension(:,:,:), intent(inout) :: psi

   double precision :: norme

   !Normalize wave function
   norme=sqrt(sum(abs(psi(:,:,:))**2)*dx*dy*dz)
   psi(:,:,:)=psi(:,:,:)/norme

   End Subroutine norm

!--- SUBROUTINE ---!

   Subroutine Proppot_3D(psi)

   Implicit none

   double complex, dimension(:,:,:), intent(inout) :: psi

   !Perform diagonal propagation
   psi=psi*Vprop

   End Subroutine Proppot_3D

!--- SUBROUTINE ---!

   Subroutine Energy3D(psi,en)

   Implicit None

   double complex, dimension(:,:,:), intent(in) :: psi

   double complex, dimension(nx,ny,nz) :: phi

   double precision, intent(out) :: en

   double precision :: ec,ep

   phi=psi !copy wave function

   !Transform to momentum representation
   call dfftw_execute_dft(planf,phi,phi)

   !Integrate < psi | k**2/2mu | psi> to get the kinetic energy
   ec=sum(rk3D(:,:,:)*((abs(phi(:,:,:)))**2))*dx*dy*dz/dble(nsize)

   !Calculate the potential part
   ep=sum(V(:,:,:)*(abs(psi(:,:,:)))**2)*dx*dy*dz

   en=ep+ec !en is the total energy (au)

   End Subroutine Energy3D

!--- SUBROUTINE ---!

   Subroutine Propkin_3D(psi)

   Implicit none

   double complex, dimension(:,:,:), intent(inout) :: psi

   !Transform to momentum representation
   call dfftw_execute_dft(planf,psi,psi)

   !Perform diagonal propagation
   psi(:,:,:)=psi(:,:,:)*ak3D(:,:,:)

   !Transform back to coordinate representation
   call dfftw_execute_dft(planb,psi,psi)

   End Subroutine Propkin_3D

!--- SUBROUTINE ---!

End Program IT3D

