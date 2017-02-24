!###########################################################################################
subroutine MakeQDDEA(nrX,ncX,ncY,XM,YM,nDMU,dmulist,dmu,dXM,dYM,RTS,qout,  &
nrD,ncDB,nrQ,ncQB,objD,AD,restD,rhsD,objQ,AQ,restQ,rhsQ,ncD,ncQ,dstart,dend)
!###########################################################################################
! Construct DDEA and QDDEA LP DUAL objects 
!########################################################################
implicit none
!input variables
integer         :: nrX,ncX,ncY
integer         :: nDMU,dmulist(nDMU),dmu
integer         :: RTS  ! 1=VRS, 2=DRS, 3=CRS, 4=IRS
integer         :: nrD,ncDB,nrQ,ncQB

real*8          :: XM(nrX,ncX),YM(nrX,ncY)
real*8          :: dXM(nrX,ncX),dYM(nrX,ncY) 
real*8          :: qout


!output variables
character*1     :: restD(nrD),restQ(nrQ)
integer         :: ncD,ncQ
integer         :: dstart,dend
real*8          :: objD(ncDB),AD(nrD,ncDB),rhsD(nrD)
real*8          :: objQ(ncQB),AQ(nrQ,ncQB),rhsQ(nrQ)

!internal variables
integer         :: ncYX
integer         :: i,j,k
integer         :: i0,i1,i2,j0,j1,j2,k1,k2
integer         :: dpick
integer         :: pointsout 
real*8          :: dx(ncX),dy(ncY)
real*8          :: IMAT(nrX,nrX)
real*8          :: plim,pmult
!####################################################################
pointsout=floor(qout*nrX)
!####################################################################
! DDEA dual model
!####################################################################
ncD=ncY+ncX

ncYX=ncY+ncX

dpick=dmulist(dmu)
dx=dXM(dmu,1:ncX)
dy=dYM(dmu,1:ncY)

objD = 0.0
objD(1:ncY) = -1.0*YM(dpick,1:ncY)
objD((ncY+1):ncYX) = XM(dpick,1:ncX)
  
restD = '<'
restD(nrX+1) = '='       ! potential super efficiency

rhsD = 0.0
rhsD(nrX+1) = 1.0

AD=0.0
AD(1:nrX,1:ncY) = YM
AD(1:nrX,(ncY+1):ncYX) = -1.0*XM

AD((nrX+1),1:ncY)=dy
AD((nrX+1),(ncY+1):ncYX)=dx
!################################ 
if(RTS== 1) then  ! 'vrs'
 objD(ncD+1)= 1.0
 AD(1:nrX,ncD+1) = -1.0 

 objD(ncD+2)= -1.0
 AD(1:nrX,ncD+2) = 1.0 

 ncD=ncD+2
end if ! on RTS=='vrs'
!################################
if(RTS== 2) then  ! 'drs'
 objD(ncD+1)= 1.0
 AD(1:nrX,ncD+1) = -1.0 

 ncD=ncD+1
end if ! on RTS=='drs'
!################################
if(RTS== 4) then  ! 'irs'
 objD(ncD+1)= -1.0
 AD(1:nrX,ncD+1) = 1.0 

 ncD=ncD+1
end if ! on RTS=='irs'
!################################

!###############################################################
! end DDEA dual section
!###############################################################


!###############################################################
! QDDEA dual model
!###############################################################
nrQ=nrX+3
ncQ=ncYX+1+nrX+1

! initialize identity matrix
IMAT=0.0
do i=1,nrX
 IMAT(i,i)=1.0
end do 
! probability level weight
plim=(1.0/float(nrX))*(float(pointsout+1)-(1E-6))
pmult=1.0/plim

objQ = 0.0
objQ(1:ncYX) = objD(1:ncYX)

restQ = '<'
restQ(nrX+1) = '=' ! potential super efficiency

rhsQ = 0.0
rhsQ(nrx+1) = 1.0

AQ = 0.0
AQ(1:(nrX+1),1:ncYX) = AD(1:(nrX+1),1:ncYX)

AQ(1:nrX,ncYX+1) = 1.0
AQ(nrX+3,ncYX+1) = -1.0

j0=ncYX+2; j2=j0+nrX-1
AQ(1:nrX,j0:j2) = -1.0*IMAT(1:nrX,1:nrX)
AQ(nrX+2,j0:j2) = 1.0/float(nrX)

j2=j2+1
AQ(nrX+2,j2) = -1.0
AQ(nrX+3,j2) = 1.0/plim
!######################################
if(RTS== 1) then  ! 'vrs'
 objQ(ncQ+1)= 1.0
 AQ(1:nrX,ncQ+1) = -1.0 

 objQ(ncQ+2)= -1.0
 AQ(1:nrX,ncQ+2) = 1.0 

 ncQ=ncQ+2
end if ! on RTS=='vrs'
!################################
if(RTS== 2) then  ! 'drs'
 objQ(ncQ+1)= 1.0
 AQ(1:nrX,ncQ+1) = -1.0 

 ncQ=ncQ+1
end if ! on RTS=='drs'
!################################
if(RTS== 4) then  ! 'irs'
 objQ(ncQ+2)= -1.0
 AQ(1:nrX,ncQ+2) = 1.0 

 ncQ=ncQ+2
end if ! on RTS=='irs'
!################################
dstart = ncYX+2
dend = dstart+nrX-1
!###############################

!################################################################
end subroutine MakeQDDEA
!################################################################


