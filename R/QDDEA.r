#############################################################################################
QDDEAR <- function(X,Y,qout=0.10,orient='in',RTS='crs',nboot=0,dmulist=NULL,DX=NULL,DY=NULL,
XList=NULL,YList=NULL,mcells=5,mlist=NULL,seedval=1001,replaceum=F,BIGM=1E5){
#############################################################################################
  replace=0
if(!is.numeric(replaceum)){
 replace=ifelse(replaceum==T,1,0)
 replace
}
seedval0=seedval
if(!is.matrix(X))  X=as.matrix(X)
if(!is.matrix(Y))  Y=as.matrix(Y)
(nrX=nrow(X));(ncX=ncol(X));(ncY=ncol(Y))
(pointsout=floor(qout*nrX))

if(orient=='ddea'& is.null(DX)|orient=='ddea'& is.null(DY)){
 stop('if orient = ddea both DX and DY must be non-null')
}


if(is.null(dmulist)) {
 dmulist=1:nrX
}
(nDMU=length(dmulist))

if(is.null(XList)) (XList=matrix(X[dmulist,],nDMU,ncX))
if(is.null(YList)) (YList=matrix(Y[dmulist,],nDMU,ncY))


if(!is.null(DX)) if(!is.matrix(DX)) DX=matrix(DX,nDMU,ncX)
if(!is.null(DY)) if(!is.matrix(DY)) DY=matrix(DY,nDMU,ncY)
 #########################
if(is.null(mlist[1])){
 # compute mlist values
 (n1=floor(1.0*sqrt(nrX)))
 (n2=floor(5.0*sqrt(nrX)))
 if(n2>0.5*nrX) n2=floor(0.5*nrX)

 mlist=floor(exp(seq(log(n1),log(n2),length.out=mcells)))
 if(max(table(mlist)>1)) mlist=floor(seq(n1,n2,length.out=mcells))

 mlist
}
(mNUM=length(mlist))
 #########################
 if(!is.null(DX)) if(ncol(DX)!=ncol(X))  stop('number of inputs  in DX do not match number of inputs in X')
 if(!is.null(DY)) if(ncol(DY)!=ncol(Y))  stop('number of outputs  in DY do not match number of outputs in Y')
 if(!is.null(DX)) if(nrow(DX)!=nDMU) stop('number of inputs  in DX do not match number of DMUs')
 if(!is.null(DY)) if(nrow(DY)!=nDMU) stop('number of outputs in DY do not match number of DMUs')

if(is.null(DX)|is.null(DY)) {
 if(orient=='in'|orient=='IN'){
  (DX=matrix(0,nDMU,ncX))
  (DY=matrix(0,nDMU,ncY))
  for(i in 1:length(dmulist)){
   DX[i,]=X[dmulist[i],]
  }
 }

  if(orient=='out'|orient=='OUT'){
   DX=matrix(0,nDMU,ncX)
   DY=matrix(0,nDMU,ncY)
   for(i in 1:length(dmulist)){
    DY[i,]=Y[dmulist[i],]
   }
  }

  if(orient=='inout'|orient=='INOUT'){
   DX=matrix(0,nDMU,ncX)
   DY=matrix(0,nDMU,ncY)
   for(i in 1:length(dmulist)){
    DX[i,]=X[dmulist[i],]
    DY[i,]=Y[dmulist[i],]
   }
  }
} # end if is.null(DX or DY)
  DX
  DY
###########################################################
(mNUM=length(mlist))
rts=3
if(RTS=='VRS'|RTS=='vrs') rts=1
if(RTS=='DRS'|RTS=='drs') rts=2
if(RTS=='IRS'|RTS=='irs') rts=4
############################################################
# create and initialize fortran output variables
  effstatus=as.integer(matrix(0,nDMU,3))
  bootstatus=as.integer(array(0,dim=c(nboot,mNUM,nDMU)))
  effvals=matrix(0,nDMU,3)
  Dprices=matrix(0,nDMU,(ncY+ncX))
  D2prices=matrix(0,nDMU,(ncY+ncX))
  boot=array(0,dim=c(nboot,mNUM,nDMU))
###########################################################

#######################################################################
# Call fortran subroutine QDDEA
#######################################################################
  Fdata =
    .Fortran('QDDEA',
             nrX=as.integer(nrX),
             ncX=as.integer(ncX),
             ncY=as.integer(ncY),
             XM=as.double(X),
             YM=as.double(Y),
             nDMU=as.integer(nDMU),
             dmulist=as.integer(dmulist),
             dXM=as.double(DX),
             dYM=as.double(DY),
             XList=as.double(XList),
             YList=as.double(YList),
             rts=as.integer(rts),
             qout=as.double(qout),
             BIGM=as.double(BIGM),
             nboot=as.integer(nboot),
             replace=as.integer(replace),
             mNUM=as.integer(mNUM),
             mlist=as.integer(mlist),
             seedval=as.integer(seedval),
             effvals=as.double(effvals),
             effstatus=as.integer(effstatus),
             Dprices=as.double(Dprices),
             D2prices=as.double(D2prices),
             boot=as.double(boot),
             bootstatus=as.integer(bootstatus),
             PACKAGE="Rqdea")

  ####################################################################

  ####################################################################
  # put Fdata results back into R objects
  effstatus=matrix(Fdata$effstatus,nDMU,3)
  bootstatus=array(as.integer(Fdata$bootstatus),dim=c(nboot,mNUM,nDMU))
  effvals=matrix(Fdata$effvals,nDMU,3)
  Dprices=matrix(Fdata$Dprices,nDMU,(ncY+ncX))
  D2prices=matrix(Fdata$D2prices,nDMU,(ncY+ncX))
  boot=array(as.integer(Fdata$boot),dim=c(nboot,mNUM,nDMU))
  seedval=Fdata$seedval
  ####################################################################
  list(effvals=effvals,effstatus=effstatus,boot=boot, bootstatus=bootstatus,
       Dprices=Dprices,D2prices=D2prices,seedval=seedval)
} # end function QDDEA
#########################################################################


