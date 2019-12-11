program virtual_packages_random
use fgsl

implicit none

real*8 :: apl,aplS,d,dS,varid,varidS,variapl,variaplS,zscored,zscoreapl,apl1,d1,apl2,sqrtapl,sqrtd,sqrtaplS,sqrtdS
integer*8 :: N,seed,i,j,nlines,nvirt,cont,nruns,nmax,nmb1,nmb2,nmbdep,nmedin,nmedout,dist
integer*8 :: i1,i2,j1,j2,cont2,kn,ko,nrew
character :: dato*128, dato2*128, dato3*128, dato4*128, salida*128, fmt1*8, flines*128
character :: virtpac*128,prpsd*128,packname1*128,packname2*128,noth*5
character, allocatable :: virtual(:,:)*128, packname(:)*128, depend(:)*128
integer*8, allocatable :: nprov(:),virtrun(:),kin(:),kout(:),packnmb(:),edgelist(:,:),kout2(:),kin2(:)
type(fgsl_rng) :: r
type(fgsl_rng_type) :: ini_fgsl

!Read from konsole the filenames. 1st is the list of virtual packages, 2nd is the list of provided by packages, 3d is the package network for
!take the node name and finally 4th is the network to analyze over runs.
call get_command_argument(1, dato)
read(dato,*)

call get_command_argument(2, dato2)
read(dato2,*)

call get_command_argument(3, dato3)
read(dato3,*)

call get_command_argument(4, dato4)
read(dato4,*)

!Label to read the number of lines
call get_command_argument(5, salida)
read(salida,*) dist
!Initialization of fgsl libraries in linux
nmax=int(1e4)
seed=21
ini_fgsl = fgsl_rng_env_setup()
fgsl_rng_default_seed=seed
ini_fgsl = fgsl_rng_default
r = fgsl_rng_alloc(ini_fgsl)


!We should read virtual packages from file
!This code tell us how many lines does the file have
call execute_command_line ('wc -l <'//trim(dato)//'>'//trim(dato)//'_'//trim(salida)//'nlines')
flines=''//trim(dato)//'_'//trim(salida)//'nlines'
open (28, file=flines, status='old') !Fichero para escribir la matriz
read(28,*) nvirt
close(28)
call execute_command_line ('rm -r '//trim(dato)//'_'//trim(salida)//'nlines')

!Once known the number of virtual packages allocate a matrix
allocate(virtual(nvirt,0:1000),nprov(nvirt),virtrun(nvirt))
nprov=0

open (28, file=dato, status='old') !File to write the matrix
do i=1, nvirt
  read(28,*) virtual(i,0)
  virtual(i,0)=trim(virtual(i,0))
enddo
close(28)

!We should read provided by packages to add to the matrix
call execute_command_line ('wc -l <'//trim(dato2)//'>'//trim(dato2)//'_'//trim(salida)//'nlines')
flines=''//trim(dato2)//'_'//trim(salida)//'nlines'
open (28, file=flines, status='old') !Fichero para escribir la matriz
read(28,*) nlines
close(28)
call execute_command_line ('rm -r '//trim(dato2)//'_'//trim(salida)//'nlines')

open (28, file=dato2, status='old') !File to write the matrix
do i=1, nlines
  read(28,*) prpsd,virtpac
  virtpac=trim(virtpac)
  prpsd=trim(prpsd)
  j=1
  do while(((virtual(j,0))/=(virtpac)).and.(j<=nvirt))
   j=j+1
  enddo
  if(j>nvirt) then
   write(*,*) 'ERROR: Packages with no virtual package associated in this distribution'
   write(*,*) trim(virtpac),' ',trim(prpsd)
   write(*,*) ' '
   write(300,*) trim(virtpac),' ',trim(prpsd)
  endif
  if(j<=nvirt) then
   nprov(j)=nprov(j)+1
   virtual(j,nprov(j))=prpsd
  endif
enddo
close(28)
do i=1, nvirt
 if(nprov(i)==0) then
 write(*,*) trim(virtual(i,0)),' ',trim(virtual(i,1)),' ',trim(virtual(i,2)),' ',trim(virtual(i,3)),' ',trim(virtual(i,4))
 endif
enddo

!Once completed the list of virtual packages
!Start average over runs
!But, first of all, we should read the list of packages and packages numbers
call execute_command_line ('wc -l <'//trim(dato3)//'>'//trim(dato3)//'_'//trim(salida)//'nlines')
flines=''//trim(dato3)//'_'//trim(salida)//'nlines'
open (28, file=flines, status='old') 
read(28,*) N
close(28)
call execute_command_line ('rm -r '//trim(dato3)//'_'//trim(salida)//'nlines')

!Once known the number of virtual packages allocate a matrix, and the degree distributions
allocate(packnmb(N),packname(N),kin(N),kout(N),edgelist(N,500),depend(100),kout2(N),kin2(N))

open (28, file=dato3, status='old') 
do i=1, N
  read(28,*) packnmb(i), packname(i)
  packname(i)=trim(packname(i))
enddo
close(28)
apl=0.0d0
aplS=0.0d0
d=0.0d0
dS=0.0d0
sqrtapl=0.0d0
sqrtaplS=0.0d0
sqrtd=0.0d0
sqrtdS=0.0d0

!Average over runs
do nruns=1, nmax
  kin=0
  kout=0
  kn=0
  edgelist=0
  !File for print the graph
  salida='graph.txt'
  open (unit = 25, file = salida)

  !We draw a random package for each virtual package I need, between 2 and nprov(i)
  do i=1,nvirt
    if (nprov(i)==1) then
      virtrun(i)=1
    else
      virtrun(i)=1+fgsl_rng_uniform_int (r,nprov(i)) !Between 1 and nprov
    endif
  enddo
 
  !Now we read the desired network, dependencies or conflicts
  call execute_command_line ('wc -l <'//trim(dato4)//'>'//trim(dato4)//'_'//trim(salida)//'nlines')
  flines=''//trim(dato4)//'_'//trim(salida)//'nlines'
  open (28, file=flines, status='old')
  read(28,*) nlines
  close(28)
  call execute_command_line ('rm -r '//trim(dato4)//'_'//trim(salida)//'nlines')
  
  open (28, file=dato4, status='old') 
  do i=1, nlines
    read(28,*) nmbdep, packname1, depend(1:nmbdep+1)
    !if (i<500) write(*,*) i, trim(packname1)
    j=1+fgsl_rng_uniform_int (r,nmbdep+1)
    packname2=depend(j)
    cont=0
    j=1
    !Try with virtual packages, for select it. Take into account that virtual packages it's only in the right hand-side of our list
    do while((cont==0).and.(j<=nvirt))
      if (packname2==virtual(j,0)) then
        packname2=virtual(j,virtrun(j))
        cont=1
      endif
      j=j+1
    enddo
    !Find for package numbers, we need it for the edgelist
    j=1
    cont=0
    do while((cont==0).and.(j<=N))
      if (packname1==packname(j)) then
        nmb1=packnmb(j)
        cont=2
      endif
      j=j+1
    enddo
    j=1
    cont=0
    do while((cont==0).and.(j<=N))
      if (packname2==packname(j)) then
        nmb2=packnmb(j)
        cont=2
      endif
      j=j+1
    enddo
    !Check if we have the link. If not, put the link 
    do j=1, kout(nmb1)
      if(nmb2==edgelist(nmb1,j)) cont=1
    enddo
    if (cont==2) then
      kin(nmb2)=kin(nmb2)+1
      kn=kn+1
      kout(nmb1)=kout(nmb1)+1
      edgelist(nmb1,kout(nmb1))=nmb2
      write(25,*) nmb1,nmb2
    endif
  enddo
  close(28)
  close(25)

  !Compute and print the Diameter and APL in each run
  call execute_command_line('Rscript R_APL > Diameter')

  open (25, file='Diameter', status='old')
  read(25,*) noth,d1
  read(25,*) noth,apl1
  
  apl=apl+apl1
  sqrtapl=sqrtapl+apl1*apl1
  d=d+d1
  sqrtd=sqrtd+d1*d1
  close(25)
  call execute_command_line('rm -r Diameter')
  call execute_command_line('rm -r graph*')

  !File for print the graph
  salida='graph.txt'
  open (unit = 25, file = salida)
  !Reshuffle network
  nrew=0
  do while(nrew<=10*kn)
    i1=1+fgsl_rng_uniform_int(r,N)
    i2=1+fgsl_rng_uniform_int(r,N)
    do while((kout(i1)==0).or.(kout(i2)==0).or.(i1==i2))
      i1=1+fgsl_rng_uniform_int(r,N)
      i2=1+fgsl_rng_uniform_int(r,N)
    enddo
    j1=1+fgsl_rng_uniform_int(r,kout(i1))
    j2=1+fgsl_rng_uniform_int(r,kout(i2))
    cont=0
    if((i1==edgelist(i2,j2)).or.(i2==edgelist(i1,j1))) cont=1
    do j=1,kout(i1)
      if (edgelist(i2,j2)==edgelist(i1,j)) cont=1
    enddo
    do j=1,kout(i2)
      if (edgelist(i1,j1)==edgelist(i2,j)) cont=1
    enddo
    if (cont==0) then
      ko=edgelist(i1,j1)
      edgelist(i1,j1)=edgelist(i2,j2)
      edgelist(i2,j2)=ko
      nrew=nrew+1
    endif
  enddo
  do i=1,N
    do j=1,kout(i)
      write(25,*) i,edgelist(i,j)
    enddo
  enddo
  close(25)

  !Compute and print the Diameter and APL in each run
  call execute_command_line('Rscript R_APL > Diameter')
  
  open (25, file='Diameter', status='old')
  read(25,*) noth,d1
  read(25,*) noth,apl1

  aplS=aplS+apl1
  sqrtaplS=sqrtaplS+apl1*apl1
  dS=dS+d1
  sqrtdS=sqrtdS+d1*d1
  close(25)
  call execute_command_line('rm -r Diameter')
  call execute_command_line('rm -r graph*')

  
  variaplS=sqrt(sqrtaplS/(1.0d0*nruns)-((aplS*aplS)/(1.0d0*nruns*nruns)))
  variapl=sqrt(sqrtapl/(1.0d0*nruns)-((apl*apl)/(1.0d0*nruns*nruns)))
  varidS=sqrt(sqrtdS/(1.0d0*nruns)-((dS*dS)/(1.0d0*nruns*nruns)))
  varid=sqrt(sqrtd/(1.0d0*nruns)-((d*d)/(1.0d0*nruns*nruns)))
  
  zscoreapl=(apl/(1.0d0*nruns))-(aplS/(1.0d0*nruns))
  zscoreapl=zscoreapl/variaplS
  zscored=(d/(1.0d0*nruns))-(dS/(1.0d0*nruns))
  zscored=zscored/varidS
  salida='APL'//trim(dato4)//'_datafile'
  open (unit = 25, file = salida)
  write(25,*) zscoreapl,(apl/(1.0d0*nruns)),(aplS/(1.0d0*nruns)),variaplS,variapl,dist,nruns  
  close(25)
  salida='Diameter'//trim(dato4)//'_datafile'
  open (unit = 25, file = salida)
  write(25,*) zscored,(d/(1.0d0*nruns)),(dS/(1.0d0*nruns)),varidS,varid,dist,nruns  
  close(25)
enddo

stop
end
