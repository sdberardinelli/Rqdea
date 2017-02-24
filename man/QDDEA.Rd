\name{QDDEA}
\alias{QDDEA}
\title{Examples for accessing QDDEA Fortran code from R}
\description{
  Call fortran QDDEA fortran code 
  }
\usage{
QDDEA(X,Y,qout=0.10,orient='in',RTS='crs',nboot=0,dmulist=NULL,DX=NULL,DY=NULL,
XList=NULL,YList=NULL,mcells=5,mlist=NULL,seedval=1001,replaceum=F,BIGM=1E5)
}

\arguments{
  \item{X}{Input matrix}
  \item{Y}{output matrix}
  \item{qout}{quantile fraction of external points}  
  \item{orient}{orientation = 'in','out','inout', 'ddea'}
  \item{RTS}{returns to scale = 'crs','vrs','drs','irs'}
  \item{nboot}{number of replications in nCm bootstrap. 0 = only return qDEA values}
  \item{dmulist}{vector of dmus to process}
  \item{DX}{matrix of input directions for dmus in dmulist}
  \item{DY}{matrix of output directions for dmus in dmulist}
  \item{XList}{matrix of input levels for dmus in dmulist}
  \item{YList}{matrix of output levels for dmus in dmulist}
  \item{mcells}{number of m levels to use in nCm bootstrap}
  \item{mlist}{user provided vector of m levels for nCm bootstrap}
  \item{seedval}{seedval}
  \item{replaceum}{sample with replacement}
  \item{BIGM}{BIGM for use in qDEA stage II}
}
\details{
  See following examples
}
\value{
  \item{effvals}{qDEA DDEA distances}
  \item{effstatus}{glpk fortran status for effvals - negative indicates not run}
  \item{boot}{array of bootstrapped efficiency scores}
  \item{bootstatus}{glpk status of values in boot - negative indicates not run}  
  \item{Dprices}{matrix of prices for DDEA model}
  \item{Dprices2}{matrix of prices for qDDEA-II model}
  \item{seedval}{seed value at completion of simulations}
}

\author{Joe Atwood}

\examples{
# not run
}
