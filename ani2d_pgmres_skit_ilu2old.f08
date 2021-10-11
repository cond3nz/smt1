program main
implicit none
real(4) t1,t2,t3,t4
integer i,n,NZ,NZN,ierr,iout,ipar(16),im,maxits,methode,verb,lendwork,leniwork,UsedWr,UsedWi,idau,idjau,idju,endidju
integer, allocatable :: ia(:),ja(:),iwork(:)
real*8 fpar(16),eps,permtol,tau1,tau2,partlur,partlurout
real*8, allocatable :: a(:),b(:),sol(:),work(:)
t1 = secnds (0.0)
methode=8
open(10,file="../sparse_matrixs/2.dat")
read(10,*) n
print *,"n =",n


!parameters for iluoo
      tau1=0.1
      tau2=sqrt(tau1)
      verb=0
      partlur=0.5

allocate(ia(n+1),b(n),sol(n))

read(10,*)(ia(i),i=1,n+1)

NZ=ia(n+1)-ia(1)
print *,"NZ = ",NZ

!len of workspaces
      lendwork=NZ*50
      leniwork=NZ*60


allocate(ja(NZ),a(NZ),work(lendwork),iwork(leniwork))

!read matrix in CSR
      read(10,*)(ja(i),i=1,NZ)
      read(10,*)(a(i),i=1,NZ)
      read(10,*)(b(i),i=1,n)

call iluoo (n,ia,ja,a,tau1,tau2,verb,work,iwork,lendwork,leniwork,partlur,partlurout,UsedWr,UsedWi,ierr)
print *," iluoo errors : ",ierr
t2 = secnds (t1)

!id massives
      
      idau=1
      idjau=n+2
      idju=idjau+n+1
      endidju=idjau+2*n

t3 = secnds (0.0)

!parameters for GMRES
      im=100
      eps  = 1.0E-08
      maxits = 10000
      iout=6
      permtol = 1.0
      ipar(2) = 1
      ipar(3) = 2
      ipar(4) = lendwork-UsedWr
      ipar(5) = im
      ipar(6) = maxits
      ipar(7)=0
      ipar(12)=0
      fpar(1) = eps
      fpar(2) = 2.22E-16

NZN=iwork(idjau+n)-iwork(idjau)
print *,"NZ precond = ",NZN

call runrc(n,b,sol,ipar,fpar,work(UsedWr:lendwork),a,ja,ia,work(idau:NZN),iwork(idjau:idjau+n),iwork(idju:endidju),methode,iout)
t4=secnds(t3)
print *, "init time = ",t2," seconds"
print *, "exec time = ",t4," seconds"
print *,"numbers of iterations = ",ipar(7)
print *,"res norm = ",fpar(6)

deallocate(ia,ja,a,b,sol,work,iwork)
end program main