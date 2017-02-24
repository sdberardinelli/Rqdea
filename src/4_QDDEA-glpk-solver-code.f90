! GLPK_subroutine_DEA-loop-DMUs-droprows
!********************************************************************
subroutine LP_GLPK(nrD,ncDB,ncD,nrQ,ncQB,ncQ,nrX,ncY,ncX,nDMU,dmulist,DYX, &
replaceAI,YXList,objD,restD,rhsD,AD,objQ,restQ,rhsQ,AQ,dstart,dend,BIGM, &
effvals,effstatus,Dprices,D2prices)
!********************************************************************
use glpk
use glpk_types
use glpk_constants
use iso_c_binding
implicit none
type(c_ptr) :: dea
type(c_ptr) :: dea2
type(c_ptr) :: qdea
!********************************************************************
!input variables
character*1 :: restD(nrD),restQ(nrQ)
integer     :: nrD,ncDB,ncD,nrQ,ncQB,ncQ
integer     :: nrX,ncY,ncX
integer     :: nDMU,dmulist(nDMU)
integer     :: replaceAI

real*8      :: DYX(nDMU,ncY+ncX)
real*8      :: YXList(nDMU,ncY+ncX)
real*8      :: objD(ncDB),rhsD(nrD),AD(nrD,ncDB)
real*8      :: objQ(ncQB),rhsQ(nrQ),AQ(nrQ,ncQB)
real*8      :: BIGM
!********************************************************************
!output variables
!DEA and QDEA outputs
integer     :: effstatus(nDMU,3)
real*8      :: effvals(nDMU,3)
real*8      :: Dprices(nDMU,ncY+ncX),D2prices(nDMU,ncY+ncX)
!********************************************************************
!internal  variables
character*3 :: objtype
character*1 :: restD2(nrD)
integer     :: ncYX,nrXp1
integer     :: irD2
integer     :: jDMU,ipick
integer     :: itmp,itmp1,itmp2,itmp3
integer     :: jtmp,jtmp1,jtmp2,jtmp3
integer     :: dstart,dend
integer     :: nchgD,nchgQ
integer     :: ijconrowD(nrD*2),ijconrowQ(nrQ*2)
integer     :: dcount,dmult(nrX)
real*8      :: obj
real*8      :: objD2(ncD)
real*8      :: rhsD2(nrD)
real*8      :: AD2(nrD,ncD)
real*8      :: solnD(ncD),solnQ(ncQ),solnD2(ncD)
real*8      :: dvals(nrX)
real*8      :: eps=10.0**(-10.0)
!********************************************************************
!glpk internal variables
integer     :: nsims

integer :: ret
integer, dimension(0:(nrD*ncD))                 :: iaD,jaD
integer, dimension(0:(nrQ*ncQ))                 :: iaQ,jaQ

integer i,j,k,ii,jj,kk,lcount
integer kD,kQ,kD2
real(kind=8), dimension(0:(nrD*ncD))            :: arD
real(kind=8), dimension(0:(nrQ*ncQ))            :: arQ

real(kind=8) :: z,z8
integer ijrowD(ncX),ijrow(ncX)
integer jcount
!********************************************************************
character(8)   :: date
character(10)  :: time
character(5)   :: zone
integer,dimension(8) :: values
real*8  :: t0,t1,t2,t3,t4,t6,t7,t8,t9,t10
real*8  :: tm0,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10
real*8  :: runtimes(20),times(20)
!********************************************************************
type(glp_smcp) :: parmQ
type(glp_smcp) :: parmD
type(glp_smcp) :: parmD2
!********************************************************************
lcount=0
objtype='min'
z8=0.0_8
!********************************************************************
call date_and_time(date,time,zone,values)
tm0=values(5)*3600+values(6)*60+values(7)+float(values(8))/1000.0

effvals=-1E5
effstatus=0

!construct the DEA data
ncYX=ncY+ncX
nrXp1=nrX+1

!turn off glpk screen messages
ret = glp_term_out(GLP_MSG_OFF)

! set up dea problem
dea  = glp_create_prob()
call glp_set_obj_dir(dea, GLP_MIN)
ret = glp_add_rows(dea, nrD)

! add rhs and row bounds
do i=1,nrD
 if(restD(i).eq.'<') then
  call glp_set_row_bnds(dea, i, GLP_UP, z8, rhsD(i))
 end if
 if(restD(i).eq.'>') then
  call glp_set_row_bnds(dea, i, GLP_LO, z8, rhsD(i))
 end if
 if(restD(i).eq.'=') then
  call glp_set_row_bnds(dea, i, GLP_FX, rhsD(i), rhsD(i))
 end if

end do

! add obj and column bounds
ret = glp_add_cols(dea, ncD)

do j=1,ncD
 call glp_set_obj_coef(dea, j, objD(j))
 call glp_set_col_bnds(dea, j, GLP_LO, z8,z8)
end do

!*************************************************
! sparce matrix construction and identification of
! change elements in dea a matrix
!*************************************************
iaD=0;jaD=0;arD=0.0
nchgD=0;ijconrowD=0;kD=0
do i=1,nrD
 if(i.ne.nrXp1) then
  do j=1,ncD
   if(abs(AD(i,j)).gt.eps) then
    kD=kD+1
    iaD(kD)=i;jaD(kD)=j;arD(kD)=AD(i,j)
!    if(i==nrD) then                         ! BIGM
!     nchgD=nchgD+1                          ! BIGM
!     ijconrowD(nchgD)=kD                    ! BIGM
!    end if                                  ! BIGM
   end if
  end do
 end if

 if(i==nrXp1) then
  do j=1,ncYX
    kD=kD+1
    iaD(kD)=i;jaD(kD)=j;arD(kD)=AD(i,j)
    nchgD=nchgD+1
    ijconrowD(nchgD)=kD
  end do
 end if

end do

call glp_load_matrix(dea, kD, iaD, jaD, arD)
call glp_init_smcp(parmD)
!solve DEA problem
!ret = glp_simplex(dea, parmD)
!obj=glp_get_obj_val(dea)
!******************************************************
! create qdea problem
qdea = glp_create_prob()

call glp_set_obj_dir(qdea, GLP_MIN)

ret = glp_add_rows(qdea, nrQ)

do i=1,nrQ
 if(restQ(i).eq.'<') then
  call glp_set_row_bnds(qdea, i, GLP_UP, z8, rhsQ(i))
 end if
 if(restQ(i).eq.'>') then
  call glp_set_row_bnds(qdea, i, GLP_LO, z8, rhsQ(i))
 end if
 if(restQ(i).eq.'=') then
  call glp_set_row_bnds(qdea, i, GLP_FX, rhsQ(i), rhsQ(i))
 end if

end do

ret = glp_add_cols(qdea, ncQ)

do j=1,ncQ
 call glp_set_obj_coef(qdea, j, objQ(j))
 call glp_set_col_bnds(qdea, j, GLP_LO, z8,z8)
end do


iaQ=0;jaQ=0;arQ=0.0
nchgQ=0;ijconrowQ=0;kQ=0
do i=1,nrQ
 if(i.ne.nrXp1) then
  do j=1,ncQ
   if(abs(AQ(i,j)).gt.eps) then
    kQ=kQ+1
    iaQ(kQ)=i;jaQ(kQ)=j;arQ(kQ)=AQ(i,j)
!    if(i==nrQ) then                          !BIGM
!     nchgQ=nchgQ+1                           !BIGM
!     ijconrowQ(nchgQ)=kQ                     !BIGM
!    end if                                   !BIGM
   end if
  end do
 end if

 if(i==nrXp1) then
  do j=1,ncYX
    kQ=kQ+1
    iaQ(kQ)=i;jaQ(kQ)=j;arQ(kQ)=AQ(i,j)
    nchgQ=nchgQ+1
    ijconrowQ(nchgQ)=kQ
  end do
 end if

end do

call glp_load_matrix(qdea, kQ, iaQ, jaQ, arQ)
call glp_init_smcp(parmQ)

!ret = glp_simplex(qdea, parmQ)
!obj=glp_get_obj_val(qdea)
! do j=1,ncQ
!  solnQ(j)=glp_get_col_prim(qdea,j)
! end do
!******************************************************

!******************************************************
! find solutions for each dmu in dmulist
jcount=0
jDMU=1
do 100 jDMU=1,nDMU
!******************************************************
!if(jDMU==34) then
! write(*,*) 'jDMU = ',jDMU
!end if

ipick=dmulist(jDMU)
!! edit objective values
!objD(1:ncYX) = -1.0*AD(ipick,1:ncYX)
!objQ(1:ncYX) = -1.0*AQ(ipick,1:ncYX)

objD(1:ncYX) = -1.0*YXList(jDMU,1:ncYX)
objQ(1:ncYX) = -1.0*YXList(jDMU,1:ncYX)

! edit directions
!AD(nrX+1,1:ncYX)=DYX(ipick,1:ncYX)
!AQ(nrX+1,1:ncYX)=DYX(ipick,1:ncYX)

AD(nrX+1,1:ncYX)=DYX(jDMU,1:ncYX)
AQ(nrX+1,1:ncYX)=DYX(jDMU,1:ncYX)

! edit unbounded restriction
!AD(nrD,1:ncYX)=-1.0*objD(1:ncYX)         ! BIGM
!AQ(nrQ,1:ncYX)=-1.0*objQ(1:ncYX)         ! BIGM


! reset objective values in glpk objects
do j=1,ncYX
 call glp_set_obj_coef(dea, j, objD(j))
end do

do j=1,ncYX
 call glp_set_obj_coef(qdea, j, objQ(j))
end do

!reset glpk dea constraint coefficients
do j=1,nchgD
 itmp=ijconrowD(j)
 arD(itmp)=AD(iaD(itmp),jaD(itmp))
end do

if(replaceAI==1) then
 arD(1:ncYX)=YXList(jDMU,1:ncYX)
end if

! reload glpk constraint coeff data
call glp_load_matrix(dea, kD, iaD, jaD, arD)

! edit qdea data

do j=1,nchgQ
 itmp=ijconrowQ(j)
 arQ(itmp)=AQ(iaQ(itmp),jaQ(itmp))
end do

if(replaceAI==1) then
 arQ(1:ncYX)=YXList(jDMU,1:ncYX)
end if

call glp_load_matrix(qdea, kQ, iaQ, jaQ, arQ)
!**************************************************

!solve DEA problem
ret = glp_simplex(dea, parmD)
ret = glp_get_status(dea)
obj=glp_get_obj_val(dea)
effvals(jDMU,1)=obj
effstatus(jDMU,1)=ret
do j=1,ncD
 solnD(j)=glp_get_col_prim(dea,j)
end do
Dprices(jDMU,1:ncYX)=solnD(1:ncYX)
!******************************************************
!QDEA stage 1
ret = glp_simplex(qdea, parmQ)
ret = glp_get_status(qdea)
effstatus(jDMU,2)=ret
if(ret .ne. 5) effstatus(jDMU,3)=ret

99 if(ret.eq.5) then
obj=glp_get_obj_val(qdea)
effvals(jDMU,2)=obj

do j=1,ncQ
  solnQ(j)=glp_get_col_prim(qdea,j)
end do

!*******************************************************
!QDEA stage 2
! find which dmu's qdea left outside the DDEA hull
dmult=0;dvals=0.0
dvals=solnQ(dstart:dend)
where(dvals>eps) dmult=1

dcount=sum(dmult)

rhsD2=rhsD
where(dmult>=1) rhsD2(1:nrX)=BIGM
!************************************************************
! set up dea2 glpk object
dea2  = glp_create_prob()

call glp_set_obj_dir(dea2, GLP_MIN)
ret = glp_add_rows(dea2, nrD)

do i=1,nrD
 if(restD(i).eq.'<') then
  call glp_set_row_bnds(dea2, i, GLP_UP, z8, rhsD2(i))
 end if
 if(restD(i).eq.'>') then
  call glp_set_row_bnds(dea2, i, GLP_LO, z8, rhsD2(i))
 end if
 if(restD(i).eq.'=') then
  call glp_set_row_bnds(dea2, i, GLP_FX, rhsD2(i), rhsD2(i))
 end if
end do

ret = glp_add_cols(dea2, ncD)

do j=1,ncD
 call glp_set_obj_coef(dea2, j, objD(j))
 call glp_set_col_bnds(dea2, j, GLP_LO, z8,z8)
end do

call glp_load_matrix(dea2, kD, iaD, jaD, arD)
call glp_init_smcp(parmD2)

ret = glp_simplex(dea2, parmD2)
ret = glp_get_status(dea2)

effstatus(jDMU,3)=ret
98 if(ret .eq. 5) then
obj=glp_get_obj_val(dea2)
effvals(jDMU,3)=obj
!write(*,*) "jDMU,objD2, statusD2 = ", jDMU,obj,ret
do j=1,ncD
 solnD2(j)=glp_get_col_prim(dea2,j)
! if(nc.lt.20) Print *, "DEA2 j,xj = ", j,solnD(j)
end do
D2prices(jDMU,1:ncYX)=solnD2(1:ncYX)

end if ! end if  98 if(ret. .eq. 5) qdea state 2
call glp_delete_prob(dea2)

!jcount=jcount+1
!if(jcount==50) then
! call cpu_time(t3)
! open(20,file='fortran-report.txt',status='unknown')
!  write(20,*) 'jDMU, runtimes = ',jDMU,times(1),times(2),times(3)
! close(20)
! jcount=0
!end if
end if   !   end if  on 99 if(ret.eq.5) then  qDEA STAGE I
100 continue

call date_and_time(date,time,zone,values)
tm10=values(5)*3600+values(6)*60+values(7)+float(values(8))/1000.0
times(4)=times(4)+tm10-tm0

!Print *, "runtime = ",runtimes(2)
! open(20,file='fortran-report.txt',status='unknown')
!  write(20,*) 'jDMU, runtimes = ',jDMU,times(1),times(2),times(3),times(4)
! close(20)

!write(*,*) 'glpk-dea-run-time =',times(1)
!write(*,*) 'glpk-qdea-run-time =',times(2)
!write(*,*) 'gplk-dea2-run-time = ',times(3)
!write(*,*) ' glpk-total-run-time =', times(1)+times(2)+times(3)
!write(*,*) 'total subroutine run time =',times(4)


!1 format(1000(1f15.6,','))
!2 format(1a1,',',1e15.6)
!3 format(2i15,1f15.6)
!4 format(1i10,',',2(1f20.8,','),1f20.8)


!open(10,file='DEBUG-QDEA-data.csv',status='unknown')
!write(10,*) 'solution found'
!write(10,1) (solnQ(i),i=1,ncQ)
!write(10,*) 'problem data'
!write(10,*) 'objective direction'
!write(10,*) 'min'
!write(10,*) 'nrows,ncolumns'
!write(10,3) nrQ,ncQ
!write(10,*) 'objective coeff'
!write(10,1) (objQ(i),i=1,ncQ)
!write(10,*) 'restQ,rhsQ'
!do i=1,nrQ
! write(10,2) restQ(i),rhsQ(i)
!end do
!write(10,*) 'Sparce matrix data'
!write(10,*) 'iaQ,jaQ,arQ'
!do i=1,kQ
! write(10,3) iaq(i),jaq(i),arQ(i)
!end do
!write(10,*) 'efficiency scores'
!do i=1,nDMU
! write(10,1) (EFFscores(i,k),k=1,3)
!end do
!close(10)


!open(10,file='EFF-scores.csv',status='unknown')
!write(10,*) 'dmu,E1,E2,E3'
!do i=1,nDMU
! write(10,4) dmulist(i),(EFFscores(i,k),k=1,3)
!end do
!close(10)

call glp_delete_prob(dea)
call glp_delete_prob(qdea)

!********************************************************************
end subroutine LP_GLPK
!********************************************************************
