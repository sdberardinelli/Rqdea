!*******************************************************************************************
subroutine QDDEA(nrX,ncX,ncY,XM,YM,nDMU,dmulist,dXM,dYM,XList,YList, &
RTS,qout,BIGM,nboot,replace,mNUM,mlist,seedval,effvals,effstatus,    &
Dprices,D2prices,boot,bootstatus)
!*******************************************************************************************
! Estimate and bootstrap DDEA and QDDEA objective values
!
!  Note: We set up DDEA and QDDEA models that allow for super efficiency
!  In doing this we allow the objective distance to be a free variable
!  with negative values indicating superefficient points.
!************************************************************************
!*************************************************************************
!input variables
integer         :: nrX,ncX,ncY
integer         :: nDMU,dmulist(nDMU)
integer         :: RTS  ! 1=VRS, 2=DRS, 3=CRS, 4=IRS
integer         :: nboot,replace
integer         :: mNUM,mlist(mNUM)
integer         :: seedval

real*8          :: XM(nrX,ncX),YM(nrX,ncY)
real*8          :: dXM(nDMU,ncX),dYM(nDMU,ncY)
real*8          :: XList(nDMU,ncX),YList(nDMU,ncY)
real*8          :: qout,BIGM
!**************************************************************************
!output variables
integer         :: effstatus(nDMU,3)
integer         :: bootstatus(nboot,mNUM,nDMU)
real*8          :: effvals(nDMU,3)
real*8          :: Dprices(nDMU,ncY+ncX),D2prices(nDMU,ncY+ncX)
real*8          :: boot(nboot,mNUM,nDMU)
!**************************************************************************
!internal variables
character*1     :: restD(nrX+1),restQ(nrX+3)
integer         :: ncDB,ncQB
integer         :: nrD,ncD,nrQ,ncQ
integer         :: dstart,dend
integer         :: dmu
integer         :: ncYX
integer         :: i,j,k,ii,jj,kk
integer         :: i0,i1,i2,j0,j1,j2,k1,k2
integer         :: replaceAI

real*8          :: dYX(nDMU,ncY+ncX)
real*8          :: objD(ncY+ncX+2),objQ(ncY+ncX+nrX+4),AD(nrX+1,ncY+ncX+2),AQ(nrx+3,ncY+ncX+nrX+4)
real*8          :: rhsD(nrX+1),rhsQ(nrx+3)
real*8          :: YXList(nDMU,ncY+ncX)
!*************************************************************************
character(8)    :: date
character(10)   :: time
character(5)    :: zone
integer,dimension(8) :: values
real*8          :: t0,t1,t2,t3,t4,t6,t7,t8,t9,t10
real*8          :: tm0,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,tm9,tm10
real*8          :: runtimes(20),times(20)
!********************************************************************
! bootstrap variables
character*1     :: restDm(:),restQm(:)
integer         :: jb,jm,nm,qpointsout
integer         :: nvals(nrX),ipick(nrX)

integer         :: nrXm,ncXm,ncYm
integer         :: nrDm,ncDBm,ncDm
integer         :: nrQm,ncQBm,ncQm
integer         :: effstatusm(:,:)

real*8          :: Xmb(:,:),Ymb(:,:)
real*8          :: objDm(:),ADm(:,:),rhsDm(:)
real*8          :: objQm(:),AQm(:,:),rhsQm(:)
real*8          :: effvalsm(:,:)



allocatable  restDm,restQm,Xmb,Ymb,objDm,ADm,rhsDm,objQm,AQm,rhsQm
allocatable  effvalsm,effstatusm
!************************************************************************
effvals = -1E5
effstatus= -1000

ncYX=ncY+ncX
nrD=nrX+1;    nrQ=nrX+3
ncDB=ncYX+2;  ncQB=ncYX+nrX+4
!*************************************************************************
dYX(1:nDMU,1:ncY)=dYM; dYX(1:nDMU,(ncY+1):(ncY+ncX))=dXM
YXList(1:nDMU,1:ncY)=YList;YXList(1:nDMU,(ncY+1):ncYX)=-1.0 * XList
!*************************************************************************
! create base lp problems
times=0.0d0

call date_and_time(date,time,zone,values)
tm1=values(5)*3600+values(6)*60+values(7)+float(values(8))/1000.0
!***************************************************************************
dmu=1
call MakeQDDEA(nrX,ncX,ncY,XM,YM,nDMU,dmulist,dmu,dXM,dYM,RTS,qout,          &
nrD,ncDB,nrQ,ncQB,objD,AD,restD,rhsD,objQ,AQ,restQ,rhsQ,ncD,ncQ,dstart,dend)
!******************
! Estimate base qDEA scores for all dmu's in dmulist
replaceAI=0

call LP_GLPK(nrD,ncDB,ncD,nrQ,ncQB,ncQ,nrX,ncY,ncX,nDMU,dmulist,DYX,         &
replaceAI,YXList,objD,restD,rhsD,AD,objQ,restQ,rhsQ,AQ,dstart,dend,BIGM,     &
effvals,effstatus,Dprices,D2prices)

!***************************************************************************


!**************************************************************************
! nCm bootstrap section
boot = -1E5
bootstatus=0

! bootstrap section
if(nboot.gt.0) then
  do i=1,nrX
    nvals(i)=i
  end do
 do jb=1,nboot
  seedval=newseed(seedval)
  call nCm(nrX,nrX,nvals,ipick,replace,seedval)

  do jm=1,mNUM  ! loop through the subsample sizes in mlist
   nm=mlist(jm)
   ! create new subsample X and Y matrix
   allocate(Xmb(nm,ncX),Ymb(nm,ncY))
   do i=1,nm
    Xmb(i,1:ncX)=XM(ipick(i),1:ncX)
    Ymb(i,1:ncY)=YM(ipick(i),1:ncY)
   end do

   !create new lp problem arrays and vectors
   nrXm=nm; ncXm=ncX; ncYm=ncY;
   nrDm=nrXm+1; ncDBm=ncYm+ncXm+2
   nrQm=nrXm+3; ncQBm=ncYm+ncXm+nrXm+4
   allocate(objDm(ncDBm),ADm(nrDm,ncDBm),restDm(nrDm),rhsDm(nrDm))
   allocate(objQm(ncQBm),AQm(nrQm,ncQBm),restQm(nrQm),rhsQm(nrQm))
   allocate(effvalsm(nDMU,3),effstatusm(nDMU,3))

   dmu=1
   call MakeQDDEA(nrXm,ncXm,ncYm,Xmb,Ymb,nDMU,dmulist,dmu,dXM,dYM,RTS,qout,         &
   nrDm,ncDBm,nrQm,ncQBm,objDm,ADm,restDm,rhsDm,objQm,AQm,restQm,rhsQm,ncDm,ncQm,   &
   dstart,dend)

   replaceAI=1
   call LP_GLPK(nrDm,ncDBm,ncDm,nrQm,ncQBm,ncQm,nrXm,ncYm,ncXm,nDMU,dmulist,DYX,    &
   replaceAI,YXList,objDm,restDm,rhsDm,ADm,objQm,restQm,rhsQm,AQm,dstart,dend,BIGM, &
   effvalsm,effstatusm,Dprices,D2prices)

   boot(jb,jm,1:nDMU)=effvalsm(1:nDMU,3)
   bootstatus(jb,jm,1:nDMU)=effstatusm(1:nDMU,3)
   deallocate(Xmb,Ymb)
   deallocate(objDm,ADm,restDm,rhsDm)
   deallocate(objQm,AQm,restQm,rhsQm)
   deallocate(effvalsm,effstatusm)
  end do ! end loop do jm=1,mNUM
 end do !end  do jb=1,nboot
end if  !  end for if(nboot.gt.0) then

!****************************************************************

end subroutine QDDEA
!****************************************************************


