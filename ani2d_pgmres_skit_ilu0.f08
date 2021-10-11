program main
implicit none
real(4) t1,t2,t3,t4
integer i,n,NZ,ierr,iout,ipar(16),m,maxits,methode
integer, allocatable :: ia(:),ja(:),jau(:),ju(:),iw(:)
real*8 fpar(16),eps,permtol
real*8, allocatable :: a(:),b(:),au(:),sol(:),wk(:)
t1 = secnds (0.0)
methode=8
open(10,file="../sparse_matrixs/4.dat")
read(10,*) n
print *,"n = ",n

!parameters for GMRES
      m   = 100
      eps  = 1.0D-9
      maxits = 10000
      iout=6
      permtol = 1.0
      ipar(2) = 1
      ipar(3) = 2
      ipar(4) = (n+3)*(m+2) + (m+1)*m/2 
      ipar(5) = m
      ipar(6) = maxits
      ipar(7)=0
      fpar(1) = eps
      fpar(2) = 1D-16


allocate(ia(n+1),b(n),iw(20*n),sol(n),wk(ipar(4)))
read(10,*)(ia(i),i=1,n+1)
NZ=ia(n+1)-ia(1)
allocate(ja(NZ),a(NZ),jau(NZ),ju(n+1),au(NZ))
read(10,*)(ja(i),i=1,NZ)
read(10,*)(a(i),i=1,NZ)
read(10,*)(b(i),i=1,n)
call ilu0 (n, a, ja, ia, au, jau, ju, iw, ierr)
t2 = secnds (t1)
print *,"ILU0 errors : ",ierr
t3 = secnds (0.0)
call runrc(n,b,sol,ipar,fpar,wk,a,ja,ia,au,jau,ju,methode,iout)
t4=secnds(t3)
print *, "init time = ",t2," seconds"
print *, "exec time = ",t4," seconds"
print *,"numbers of iterations = ",ipar(7)
print *,"residual/error norm = ",fpar(6)
deallocate(ia,ja,a,b,jau,ju,iw,au,sol,wk)
end program main
