!*********************************************************
! Miller's Fortran90 Random Number Generation Code
!*********************************************************
! A Fortran-program for MT19937: Real number version
!
! Code converted using TO_F90 by Alan Miller
! Date: 1999-11-26  Time: 17:09:23
! Latest revision - 5 February 2002
! A new seed initialization routine has been added based upon the new
! C version dated 26 January 2002.
! This version assumes that integer overflows do NOT cause crashes.
! This version is compatible with Lahey's ELF90 compiler,
! and should be compatible with most full Fortran 90 or 95 compilers.
! Notice the strange way in which umask is specified for ELF90.
 
! genrand() generates one pseudorandom real number (double) which is
! uniformly distributed on [0,1]-interval, for each call.
! sgenrand(seed) set initial values to the working area of 624 words.
! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
! integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.

! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Library General Public License as published by
! the Free Software Foundation; either version 2 of the License, or (at your
! option) any later version.   This library is distributed in the hope that
! it will be useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this library; if not, write to the Free Foundation, Inc.,
! 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.

!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed

! This program uses the following standard intrinsics.
!   ishft(i,n): If n > 0, shifts bits in i by n positions to left.
!               If n < 0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

!***********************************************************************

!***********************************************************************

MODULE rgen
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! Period parameters
INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681
!                                    constant vector a
INTEGER, PARAMETER :: umask = -2147483647 - 1
!                                    most significant w-r bits
INTEGER, PARAMETER :: lmask =  2147483647
!                                    least significant r bits
! Tempering parameters
INTEGER, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

!                     the array for the state vector
INTEGER, SAVE      :: mt(0:n-1), mti = n1
!                     mti==N+1 means mt[N] is not initialized

PRIVATE
PUBLIC :: dp, sgrnd, grnd, init_genrand
!PUBLIC :: ugen,normgen

CONTAINS


SUBROUTINE sgrnd(seed)
! This is the original version of the seeding routine.
! It was replaced in the Japanese version in C on 26 January 2002
! It is recommended that routine init_genrand is used instead.

INTEGER, INTENT(IN)   :: seed

!    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

mt(0)= IAND(seed, -1)
DO  mti=1,n-1
  mt(mti) = IAND(69069 * mt(mti-1), -1)
END DO

RETURN
END SUBROUTINE sgrnd
!***********************************************************************

SUBROUTINE init_genrand(seed)
! This initialization is based upon the multiplier given on p.106 of the
! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.

! This version assumes that integer overflow does NOT cause a crash.

INTEGER, INTENT(IN)  :: seed

INTEGER  :: latest

mt(0) = seed
latest = seed
DO mti = 1, n-1
  latest = IEOR( latest, ISHFT( latest, -30 ) )
  latest = latest * 1812433253 + mti
  mt(mti) = latest
END DO

RETURN
END SUBROUTINE init_genrand
!***********************************************************************


!************************************************************************
FUNCTION grnd() RESULT(fn_val)
real*8    :: fn_val

INTEGER, SAVE :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * MATA for x=0,1
INTEGER       :: kk, y

! These statement functions have been replaced with separate functions
! tshftu(y) = ISHFT(y,-11)
! tshfts(y) = ISHFT(y,7)
! tshftt(y) = ISHFT(y,15)
! tshftl(y) = ISHFT(y,-18)

IF(mti >= n) THEN
!                       generate N words at one time
  IF(mti == n+1) THEN
!                            if sgrnd() has not been called,
    CALL sgrnd(4357)
!                              a default initial seed is used
  END IF
  
  DO  kk = 0, n-m-1
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  DO  kk = n-m, n-2
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  y = IOR(IAND(mt(n-1),umask), IAND(mt(0),lmask))
  mt(n-1) = IEOR(IEOR(mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
  mti = 0
END IF

y = mt(mti)
mti = mti + 1
y = IEOR(y, tshftu(y))
y = IEOR(y, IAND(tshfts(y),tmaskb))
y = IEOR(y, IAND(tshftt(y),tmaskc))
y = IEOR(y, tshftl(y))

IF(y < 0) THEN
  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
ELSE
  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
END IF

RETURN
END FUNCTION grnd


FUNCTION tshftu(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-11)
RETURN
END FUNCTION tshftu


FUNCTION tshfts(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,7)
RETURN
END FUNCTION tshfts


FUNCTION tshftt(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,15)
RETURN
END FUNCTION tshftt


FUNCTION tshftl(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-18)
RETURN
END FUNCTION tshftl
!***************************************************************************

!************************************************************
END MODULE rgen

!*************************************************************


!***************************************************************************
subroutine ugen(n,seedval,x)
!written/modified by Joe Atwood
! 
 use rgen   
 implicit none
 integer n,seedval
 integer i 
 real*8    x(n) 
 CALL init_genrand(seedval)
 do i=1,n
  x(i)=grnd()
 end do
end subroutine ugen
!**************************************************************

!****************************************************************************!      
subroutine normgen(nobs,seedum,E)                                           !
!written/modified by Joe Atwood
 implicit none                                                         !
 ! Variables                                                            !
 INTEGER nobs,seedum                                             !
 real*8 E(nobs)
 integer seed
 real*8 uniformRand1(nobs), uniformRand2(nobs)                          !
 real*8 pi                                                        !
 !! **** NORMAL  DRAWS **** !                                                 !
 if(seedum.eq.0) seedum=1
 seed=seedum                                                            !
 call ugen(nobs,seed,uniformRand1)                                       !
 seed=-1*seedum
 call ugen(nobs,seed,uniformRand2)                                       !
 pi = 4.0d0 * atan ( 1.0d0 )                                            !
 E=sqrt(-2.0d0*log(uniformRand1))*cos(2.0d0*pi*uniformRand2)            !
end subroutine normgen                                                   !
 !****************************************************************************!   

!***********************************************************
function newseed(seedval)
! written by Joe Atwood
 implicit none
 integer seedval,newseed
 real*8  u(2)

 call ugen(2,seedval,u)
 newseed=floor(u(1)*2147483645)
 if(newseed==0) newseed=floor(u(2)*2147483645) 
 if(newseed==0) newseed=1
end function newseed
!***********************************************************

!**********************************************
subroutine nCm(n,m,nvals,mvals,replace,seedval)
!**********************************************
! written by Joe Atwood
 implicit none
 !input variables
 integer n,m,nvals(n),replace,seedval
 !output variables
 integer mvals(m)
 !internal variables
 integer ivals(n),nleft,i,k,pick(m)
 real*8 u(m)
!***********
 ivals=nvals
 call ugen(m,seedval,u)

 if(replace>0) then
  pick=floor(u*float(n))+1
  where(pick>n) pick=n
  mvals=ivals(pick)
 end if ! on if(replace==0)

 if(replace==0) then
  do i=1,m
   nleft=n-i+1
   k=floor(u(i)*float(nleft))+1
   if(k>nleft) k=nleft
   mvals(i)=ivals(k)
   ivals(k)=ivals(nleft)
  end do
 end if ! on if(replace>0)

!************************************************************************
end subroutine nCm
!************************************************************************
