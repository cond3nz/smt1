program main
implicit none
real(4) t1,t2,t3,t4
integer i,n,NZ,NZN,ierr,iout,ipar(16),im,maxits,methode,verb,lendwork,leniwork,UsedWr,UsedWi,luv,mi,ilu,iu,ir
integer lui,u,mv,v
integer, allocatable :: ia(:),ja(:),iwork(:)
real*8 fpar(16),eps,tau1,tau2,partlur,partlurout
real*8, allocatable :: a(:),b(:),sol(:),work(:),y(:)
t1 = secnds (0.0)
methode=8
open(10,file="../sparse_matrixs/3.dat") 
read(10,*) n
print *,"n =",n


!parameters for iluoo
      tau1=0.1
      tau2=10*sqrt(tau1)
      verb=1
      partlur=0.5

allocate(ia(n+1),b(n),sol(n),y(n))

read(10,*)(ia(i),i=1,n+1)

NZ=ia(n+1)-ia(1)
print *,"NZ = ",NZ

!len of workspaces
      lendwork=NZ*500
      leniwork=NZ*600


allocate(ja(NZ),a(NZ),work(lendwork),iwork(leniwork))

!read matrix in CSR
      read(10,*)(ja(i),i=1,NZ)
      read(10,*)(a(i),i=1,NZ)
      read(10,*)(b(i),i=1,n)

call iluoo (n,ia,ja,a,tau1,tau2,verb,work,iwork,lendwork,leniwork,partlur,partlurout,UsedWr,UsedWi,ierr)
print *," iluoo errors : ",ierr
t2 = secnds (t1)


t3 = secnds (0.0)

!parameters for GMRES
      im=100
      eps  = 1E-10
      maxits = 10000
      ipar(2) = 1
      ipar(3) = 2
      ipar(4) = lendwork-UsedWr
      ipar(5) = im
      ipar(6) = maxits
      ipar(7)=0
      fpar(1) = eps
      fpar(2) = 1E-18

!id massives
      u    = 1
      v    = u    + n
      mv   = v    + n 
      luv  = mv   + n

      mi   = 1
      ilu  = mi   + n+1
      iu   = ilu  + n+1
      ir   = iu   + n
      lui  = ir   + n 

      NZN=UsedWr-luv

print *,"NZ precond = ",NZN

!call runrc(n,b,sol,ipar,fpar,work(UsedWr:lendwork),a,ja,ia,work(luv:luv+NZN-1),iwork(lui:lui+NZN-1),iwork(ilu:ir-1),methode)
t4=secnds(t3)
print *, "init time = ",t2," seconds"
print *, "exec time = ",t4," seconds"
print *,"numbers of iterations = ",ipar(7)
print *,"error code of MGSRO = ",ipar(12)
print *,"number of initializations = ",ipar(13)
print *,"residual norm = ",fpar(6)



deallocate(ia,ja,a,b,sol,work,iwork)
end program main