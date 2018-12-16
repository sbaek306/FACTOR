! simulation for Suficient dimention reduction with censoring (SDRC) project
! following the procedure provided in paper SDRC

! functions
module SDRCFunctionsVar
	! definition of variables
	! parameter define:
	double precision :: exp1, randval, epsilonval
  double precision :: tempconsb=0.05,tempconshd=0.5,tempconshn=0.02
  double precision :: tempconsed=0.05,tempconsen=0.02,tempconshn2=0.15
  integer, parameter :: tt=1000
	double precision :: mean= 0.0, r1 = 0.5, kval=1.0, std= 1.0,pi=3.1415629, bandb, bandh, bandhe, censors=450
	! double precision, dimension(10), allocatable :: outpcts
	! character (len=10) :: programtype, curve, error
	 character (len=3) :: filenum1,filenum2,filenum3,filenum4
	! filename argument
	integer, parameter :: out_unit1 = 1, out_unit2 = 2, out_unit3 = 3, out_unit4 = 4
	integer, parameter :: samplen=300, d=1, p=5, p1=d*(p-d),lw=p1*(p1*3+13),lwt = 2*(2*3+13)
	double precision, dimension(samplen,p):: xall
	double precision, dimension(samplen,samplen):: zdiffall
	double precision, dimension(samplen) :: tall, logt, callval, zall, deltainput
	double precision, dimension(:,:), allocatable:: parasTrue, paras0
	double precision, dimension(d*(p-d),d*(p-d)) :: score2Val, score2ValInv
	character (len=4) :: typeKer = "quar"
	character (len=5) :: typeKer1 = "quar1",typeker2 = "quar2"
	character(len=8) :: fmtname
end module SDRCFunctionsVar
! generate simulation data
! generate required distributions


! Loading data 
subroutine getx(iter,p,n,x)
  integer,intent(in) :: iter,p,n
  double precision,dimension(iter*p,n),intent(out) :: x
  integer(2) :: i,j
  open(11,file='F.dat')
  do i=1,iter*p
     read(11,*) (x(i,j),j=1,n)
  end do
  close(11)
  return
end subroutine getx

subroutine getz(iter,n,z)
  integer,intent(in) :: iter,n
  double precision,dimension(iter,n),intent(out) :: z
  open(11,file='Ycc.dat')
  do i=1,iter
     read(11,*) (z(i,j),j=1,n)
  end do
  close(11)
  return
end subroutine getz

subroutine getd(iter,n,d)
  integer,intent(in) :: iter,n
  double precision,dimension(iter,n),intent(out) :: d
  open(11,file='Dcc.dat')
  do i=1,iter
     read(11,*) (d(i,j),j=1,n)
  end do
  close(11)
  return
end subroutine getd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Main program !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program sdrc
! the main program for sdrc
use SDRCFunctionsVar
implicit none
	integer :: eflag, i, j, t, e, ii
  integer :: a1,a2
	double precision, dimension(2) :: tempx, fvalx
	double precision :: temp, ntemp, xstd, c, criterion, vic
	double precision :: tol = 1e-8
	double precision,dimension(p1) :: fval
	double precision,dimension(lw) :: wv
	double precision,dimension(lwt) :: wvt
	double precision, dimension(d*(p-d)) :: scoreVal,estVar1
	double precision, dimension(d*p) :: parasInitial
	double precision, dimension((d+1)*(p-d-1)) :: scoreValnew1,scoreValnew2
    double precision, dimension(d*(p-d),d*(p-d)) :: estVarMat
	external  :: scoreOne, testfun
	double precision,dimension(samplen) :: xvar, randx, callval1
	double precision, dimension(tt*p,samplen) :: bigx
	double precision, dimension(tt,samplen) :: bigz
	double precision, dimension(tt,samplen) :: bigd
  double precision, dimension(tt,p-d,d) :: est 
  double precision, dimension(tt,p-d) :: Var
  double precision, dimension(tt,p-d,p-d) :: VarMat

	call getx(tt,p,samplen,bigx)
	call getz(tt,samplen,bigz)
	call getd(tt,samplen,bigd)

  allocate(parasTrue((p-d),d), paras0((p-d),d))
  parasTrue(1,1)= 0.0
  parasTrue(2,1)= 0.3 
  parasTrue(3,1)= -0.3 
  parasTrue(4,1)= -1.0
!  parasTrue(5,1)= -0.3 
!  parasTrue(6,1)= 0.0 
!  parasTrue(7,1)= -1.0 
!  parasTrue(8,1)= 0.5
!  parasTrue(9,1)= -0.5 
!  parasTrue(10,1)= 0.0 

  open(1,file='d31cc_est.txt')
  open(2,file='d31cc_var.txt')
!  open(3,file='a8_varMat.txt')
  
  do ii = 1,tt
      a1=(ii-1)*p+1
      a2=ii*p
      xall=transpose(bigx(a1:a2,:)) ! n by p                                                      
      zall=bigz(ii,:)
      do i = 1,samplen
         do j = 1,samplen
            zdiffall(i,j) = zall(j)-zall(i)
         enddo
      enddo
      deltainput=bigd(ii,:)
    
!    print*,xall(5,1),zdiffall(5,1),zall(5),deltainput(5)
	 ntemp = samplen**(-1./5.)
	 xvar = 0.0
	 do i = 1,samplen
        xvar(i) = xall(i,d)
        do j = 1,(p-d)
          xvar(i) = xvar(i)+xall(i,d+j)*parasTrue(j+1,d)
        enddo
   enddo
   xstd = sqrt(sum((xvar-sum(xvar)/dble(samplen))**2)/dble(samplen))
   bandh = xstd*ntemp
   bandhe = xstd*ntemp
   bandb = xstd*ntemp

!call init_random_seed() ! by SC
    do i = 1,p-d
      do j = 1,d
!call RANDOM_NUMBER(epsilonval) ! by SC
        paras0(i,j) = parasTrue(i,j) + (rand(0)-0.5)/4.!(epsilonval-0.5)*2 ! by SC
      enddo
    enddo

  	print *, "Working on",ii,"th simulation....."
	 call scoreOne(d*(p-d),paras0,scoreVal,eflag)
	 call hybrd10(scoreOne,d*(p-d),paras0,fval,tol,eflag,wv,lw)
   call scoreTwo(d*(p-d),paras0,estVar1,estVarMat)

      est(ii,:,:)=paras0
      Var(ii,:)=estVar1
!      VarMat(ii,:,:)=estVarMat

	    write (1,*) est(ii,:,:)
      write (2,*) Var(ii,:)
!      write (3,*) VarMat(ii,:,:)

	enddo

  close(1)
  close(2)
!  close(3)

end program sdrc

subroutine scoreOne(n,paras1,scoreVal,eflag)
! paras: parameter \beta in equation (5)
! xinput, deltainput, zinput: matrix X, \Delat vector, Z vector in equation (5)
! Y is indicator function which calculated in this function
! samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1

  use SDRCFunctionsVar
  integer :: i, j, k, n, t
  double precision, dimension(d):: lambdaNum,temp2,termOne,termTwo,lambdaVal
  double precision, dimension(p-d):: expectNum
  double precision, dimension(samplen,samplen,d):: btxall
  double precision, dimension(samplen,d)::lambdaVec,btransx,khprimevec,numeVec
  double precision, dimension(samplen,(p-d)):: expectVec
  double precision:: lambdaDen, expectDen, temp1, eps 
  double precision, dimension((p-d),d):: paras1
  double precision, dimension(n), intent(out):: scoreVal
  integer,intent(in) :: eflag
  double precision, dimension(samplen)::  betaXtemp2,  denoVec, val2!vector of denominator in hazard estimation

  eps = 1e-8
  scoreVal = 0.0
  do i = 1,samplen
    do j = 1,d
      btransx(i,j) = xall(i,j)
      do t = 1,(p-d)
        btransx(i,j) = btransx(i,j)+paras1(t,j)*xall(i,d+t)
      enddo
    enddo
  enddo
  do i = 1,samplen
    do j = 1,samplen
      btxall(i,j,:) = btransx(j,:)-btransx(i,:)
    enddo
  enddo
  do i = 1,samplen
    call kernel(btxall(i,:,:)/bandh/tempconshd,samplen,betaXtemp2,'norm',d)
    betaXtemp2 = betaXtemp2/bandh/tempconshd !kernel values
    do t = 1,samplen
      temp1 = 0.0
      do j = 1,samplen
        if (zall(j)>=zall(t)) then
          temp1 = temp1 + betaXtemp2(j)
        end if
      enddo
      denoVec(t) = temp1 + eps ! by SC
    enddo
    call kernel(zdiffall(i,:)/bandb/tempconsb,samplen,val2,'norm',1)
    val2 = val2/bandb/tempconsb
    lambdaDen =  sum(val2*deltainput*betaXtemp2/denoVec)
!print*,lambdaDen
    call kernel2(btxall(i,:,:)/bandh/tempconshn,samplen,khprimevec)!kernel values of K\prime(.)
    khprimevec = khprimevec/(bandh**2*tempconshn**2)
    do t = 1,samplen
      temp2 = 0.0
      do j = 1,samplen
        if (zall(j)>=zall(t)) then
          temp2 = temp2 + khprimevec(j,:)
        end if
      enddo
      numeVec(t,:) = temp2
    enddo
    do t = 1,d
      termOne(t) = sum(val2*deltainput*khprimevec(:,t)/denoVec)
      termTwo(t) = sum(val2*deltainput*betaXtemp2*numeVec(:,t)/(denoVec*denoVec))
      lambdaNum(t) = -termOne(t)+termTwo(t)
    enddo

    do t = 1,d
      lambdaVec(i,t) = lambdaNum(t)/(lambdaDen + eps) ! by SC
    enddo
! print*,"lambda",lambdavec 
    
    call kernel(btxall(i,:,:)/bandhe/tempconsed,samplen,val2,'norm',d)
    do j = 1,samplen
      if (zdiffall(i,j)>=0) then
        betaXtemp2(j) = val2(j)
      else
        betaXtemp2(j) = 0.0
      end if
    enddo
    expectDen = sum(betaXtemp2)/sum(val2)
        call kernel(btxall(i,:,:)/bandhe/tempconsen,samplen,val2,'norm',d)
        do j = 1,samplen
            if (zdiffall(i,j)>=0) then
                betaXtemp2(j) = val2(j)
            else
                betaXtemp2(j) = 0.0
            end if
        enddo

    do t = 1,(p-d)
      expectNum(t) = sum(betaXtemp2*xall(:,d+t))/(sum(val2) + eps) ! by SC
    enddo

    do j = 1,(p-d)
      expectVec(i,j) = xall(i,d+j)-expectNum(j)/(expectDen + eps) ! by SC
    enddo
  enddo
! print*,"expectation",expectVec
  do i = 1,d
    do j = 1,(p-d)
      do k = 1,samplen
        scoreVal((i-1)*(p-d)+j) = scoreVal((i-1)*(p-d)+j)+&
                    deltainput(k)*lambdaVec(k,i)*expectVec(k,j)
      enddo
    enddo
  enddo
!  print*,deltainput
!  print*,lambdaVec
!  print*,expectNum
!  print*,expectVec
!  print*,xall
!  print*,scoreVal
  scoreVal = scoreVal/samplen
!  print*,scoreVal

end subroutine scoreOne


subroutine scoreTwo(n,paras1,estimateVar,scoreValMatInv)
! calculate the inference part: dM=dN+Y(s), among which the Y(s) part.
! paras: parameter \beta in equation (5)
! xinput, deltainput, zinput: matrix X, \Delat vector, Z vector in equation (5)
! Y is indicator function which calculated in this function
! samplen: number of observations, p: dimension of all parameters, d: dimension of upper block, where beta's are 1
  use SDRCFunctionsVar
  integer :: i, j, k, n, t, ii,jj
  double precision, dimension(n,n), intent(out):: scoreValMatInv
  double precision, dimension(n), intent(out):: estimateVar
  double precision, dimension(n,n):: scoreValMat
  double precision, dimension(samplen,n,n):: scoreValMatall
  double precision, dimension(samplen,n):: meanVec
  double precision, dimension(d):: lambdaNum, temp2, termOne, termTwo,lambdaVal
  double precision, dimension(p-d):: expectNum
  double precision, dimension(samplen,samplen,d):: btxall
  double precision, dimension(samplen,d):: lambdaVec, btransx,khprimevec,numeVec
  double precision, dimension(samplen,(p-d)):: expectVec
  double precision:: lambdaDen, expectDen, temp1,eps
  double precision, dimension((p-d),d), intent(in):: paras1
  double precision, dimension(samplen)::  betaXtemp2,  denoVec, val2!vector of denominator in hazard estimation
  eps=1e-8
  scoreValMat = 0.0
  do i = 1,samplen
    do j = 1,d
      btransx(i,j) = xall(i,j)
      do t = 1,(p-d)
        btransx(i,j) = btransx(i,j)+paras1(t,j)*xall(i,d+t)
      enddo
    enddo
  enddo
  do i = 1,samplen
    do j = 1,samplen
      btxall(i,j,:) = btransx(j,:)-btransx(i,:)
    enddo
  enddo
  do i = 1,samplen
    call kernel(btxall(i,:,:)/bandh/tempconshd,samplen,betaXtemp2,'norm',d)
    betaXtemp2 = betaXtemp2/bandh/tempconshd !kernel values
    do t = 1,samplen
      temp1 = 0.0
      do j = 1,samplen
        if (zall(j)>=zall(t)) then
          temp1 = temp1 + betaXtemp2(j)
        end if
      enddo
      denoVec(t) = temp1 + eps ! +eps...by SC
    enddo
!print*,betaXtemp2
    call kernel(btxall(i,:,:)/bandh/tempconshn,samplen,betaXtemp2,'norm',d)
!    print*,bandh,bandb,bandhe
!    print*,btxall(:,:,1)
!    print*,betaXtemp2
    betaXtemp2 = betaXtemp2/bandh/tempconshn !kernel values
    call kernel(zdiffall(i,:)/bandb/tempconsb,samplen,val2,'norm',1)
!    print*,val2
    val2 = val2/bandb/tempconsb
    lambdaDen =  sum(val2*deltainput*betaXtemp2/denoVec)
    call kernel2(btxall(i,:,:)/bandh/tempconshn2,samplen,khprimevec)!kernel values of K\prime(.)
!    print*,khprimevec
    khprimevec = khprimevec/(bandh*bandh*tempconshn2*tempconshn2)
    do t = 1,samplen
      temp2 = 0.0
      do j = 1,samplen
        if (zall(j)>=zall(t)) then
          temp2 = temp2 + khprimevec(j,:)
        end if
      enddo
      numeVec(t,:) = temp2
    enddo
    do t = 1,d
      termOne(t) = sum(val2*deltainput*khprimevec(:,t)/denoVec)
      termTwo(t) = sum(val2*deltainput*betaXtemp2*numeVec(:,t)/(denoVec*denoVec))
      lambdaNum(t) = -termOne(t)+termTwo(t)
    enddo

!print*,"val2 sum",sum(val2) !okay
!print*,"deltainput sum",deltainput(1) !okay
!print*,"khprimevec",khprimevec
!print*,"betaXtemp2 sum",betaXtemp2
!print*,"numeVec",numeVec(1,1)
!print*,"denoVec",denoVec!okay

    do t = 1,d
      lambdaVec(i,t) = lambdaNum(t)/(lambdaDen + eps) ! by SC
    enddo
! print*,"lambda",lambdaNum
  
    call kernel(btxall(i,:,:)/bandhe/tempconsed,samplen,val2,'norm',d)
    val2 = val2/bandhe/tempconsed
    do j = 1,samplen
      if (zdiffall(i,j)>=0) then
        betaXtemp2(j) = val2(j)
      else
        betaXtemp2(j) = 0.0
      endif
    enddo
    expectDen = sum(betaXtemp2)
    call kernel(btxall(i,:,:)/bandhe/tempconsen,samplen,val2,'norm',d)
    val2 = val2/bandhe/tempconsen
    do j = 1,samplen
      if (zdiffall(i,j)>=0) then
        betaXtemp2(j) = val2(j)
      else
        betaXtemp2(j) = 0.0
      endif
    enddo
    do t = 1,(p-d)
      expectNum(t) = sum(betaXtemp2*xall(:,d+t))
    enddo

    do j = 1,(p-d)
      expectVec(i,j) = xall(i,d+j)-expectNum(j)/(expectDen + eps) ! by SC
    enddo
  enddo
! print*,"expectation",expectVec
  do i = 1,d
    do j = 1,(p-d)
      do k = 1,samplen
        meanVec(k,(i-1)*(p-d)+j) = lambdaVec(k,i)*expectVec(k,j)
      enddo
    enddo
  enddo
  do i = 1,n
    do j = 1,n
      do k = 1,samplen
        scoreValMatall(k,i,j) = deltainput(k)*meanVec(k,i)*meanVec(k,j)
      enddo
    enddo
  enddo
  do i = 1,n
    do j = 1,n
      scoreValMat(i,j) = sum(scoreValMatall(:,i,j))
    enddo
  enddo
  scoreValMat = scoreValMat/samplen
!  do i = i,n  ! by SC
!    scoreValMat(i,i)=scoreValMat(i,i) + eps !by SC
!  end do ! by SC
  call matInv(scoreValMat,scoreValMatInv,n)
  do i = 1,n
    estimateVar(i) = scoreValMatInv(i,i)
  enddo

! print*,"mat:",estimateVar
  !print*,tempconsh,tempconsb,tempconse
  !print*,"inverse:",scoreValMatInv
end subroutine scoreTwo

!!!!!!!!!!!!!!!!!!!!!!!
!!!!! subroutines !!!!!
!!!!!!!!!!!!!!!!!!!!!!!

subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
      integer n,ldfjac,iflag,ml,mu
      double precision epsfcn
      double precision x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)

      integer i,j,k,msum
      double precision eps,epsmch,h,temp,zero
      double precision dpmpar
      data zero /0.0d0/
  external fcn
      epsmch = dpmpar(1)

      eps = dsqrt(dmax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum .lt. n) go to 40

         do 20 j = 1, n
            temp = x(j)
            h = eps*dabs(temp)
            if (h .eq. zero) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag)
            if (iflag .lt. 0) go to 30
            x(j) = temp
            do 10 i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
   10          continue
   20       continue
   30    continue
         go to 110
   40 continue
         do 90 k = 1, msum
            do 60 j = k, n, msum
               wa2(j) = x(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               x(j) = wa2(j) + h
   60          continue
            call fcn(n,x,wa1,iflag)
            if (iflag .lt. 0) go to 100
            do 80 j = k, n, msum
               x(j) = wa2(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               do 70 i = 1, n
                  fjac(i,j) = zero
                  if (i .ge. j - mu .and. i .le. j + ml)  fjac(i,j) = (wa1(i) - fvec(i))/h
   70             continue
   80          continue
   90       continue
  100    continue
  110 continue
      return

end subroutine

subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,&
    mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,&
    qtf,wa1,wa2,wa3,wa4)
  use SDRCFunctionsVar
  integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
  double precision xtol,epsfcn,factor
  double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),&
    qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
  external fcn
!c     **********
!c
!c     subroutine hybrd
!c
!c     the purpose of hybrd is to find a zero of a system of
!c     n nonlinear functions in n variables by a modification
!c     of the powell hybrid method. the user must provide a
!c     subroutine which calculates the functions. the jacobian is
!c     then calculated by a forward-difference approximation.
!c
!c     the subroutine statement is
!c
!c       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
!c                        diag,mode,factor,nprint,info,nfev,fjac,
!c                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
!c
!c     where
!c
!c       fcn is the name of the user-supplied subroutine which
!c         calculates the functions. fcn must be declared
!c         in an external statement in the user calling
!c         program, and should be written as follows.
!c
!c         subroutine fcn(n,x,fvec,iflag)
!c         integer n,iflag
!c         double precision x(n),fvec(n)
!c         ----------
!c         calculate the functions at x and
!c         return this vector in fvec.
!c         ---------
!c         return
!c         end
!c
!c         the value of iflag should not be changed by fcn unless
!c         the user wants to terminate execution of hybrd.
!c         in this case set iflag to a negative integer.
!c
!c       n is a positive integer input variable set to the number
!c         of functions and variables.
!c
!c       x is an array of length n. on input x must contain
!c         an initial estimate of the solution vector. on output x
!c         contains the final estimate of the solution vector.
!c
!c       fvec is an output array of length n which contains
!c         the functions evaluated at the output x.
!c
!c       xtol is a nonnegative input variable. termination
!c         occurs when the relative error between two consecutive
!c         iterates is at most xtol.
!c
!c       maxfev is a positive integer input variable. termination
!c         occurs when the number of calls to fcn is at least maxfev
!c         by the end of an iteration.
!c
!c       ml is a nonnegative integer input variable which specifies
!c         the number of subdiagonals within the band of the
!c         jacobian matrix. if the jacobian is not banded, set
!c         ml to at least n - 1.
!c
!c       mu is a nonnegative integer input variable which specifies
!c         the number of superdiagonals within the band of the
!c         jacobian matrix. if the jacobian is not banded, set
!c         mu to at least n - 1.
!c
!c       epsfcn is an input variable used in determining a suitable
!c         step length for the forward-difference approximation. this
!c         approximation assumes that the relative errors in the
!c         functions are of the order of epsfcn. if epsfcn is less
!c         than the machine precision, it is assumed that the relative
!c         errors in the functions are of the order of the machine
!c         precision.
!c
!c       diag is an array of length n. if mode = 1 (see
!c         below), diag is internally set. if mode = 2, diag
!c         must contain positive entries that serve as
!c         multiplicative scale factors for the variables.
!c
!c       mode is an integer input variable. if mode = 1, the
!c         variables will be scaled internally. if mode = 2,
!c         the scaling is specified by the input diag. other
!c         values of mode are equivalent to mode = 1.
!c
!c       factor is a positive input variable used in determining the
!c         initial step bound. this bound is set to the product of
!c         factor and the euclidean norm of diag*x if nonzero, or else
!c         to factor itself. in most cases factor should lie in the
!c         interval (.1,100.). 100. is a generally recommended value.
!c
!c       nprint is an integer input variable that enables controlled
!c         printing of iterates if it is positive. in this case,
!c         fcn is called with iflag = 0 at the beginning of the first
!c         iteration and every nprint iterations thereafter and
!c         immediately prior to return, with x and fvec available
!c         for printing. if nprint is not positive, no special calls
!c         of fcn with iflag = 0 are made.
!c
!c       info is an integer output variable. if the user has
!c         terminated execution, info is set to the (negative)
!c         value of iflag. see description of fcn. otherwise,
!c         info is set as follows.
!c
!c         info = 0   improper input parameters.
!c
!c         info = 1   relative error between two consecutive iterates
!c                    is at most xtol.
!c
!c         info = 2   number of calls to fcn has reached or exceeded
!c                    maxfev.
!c
!c         info = 3   xtol is too small. no further improvement in
!c                    the approximate solution x is possible.
!c
!c         info = 4   iteration is not making good progress, as
!c                    measured by the improvement from the last
!c                    five jacobian evaluations.
!c
!c         info = 5   iteration is not making good progress, as
!c                    measured by the improvement from the last
!c                    twenty iterations.
!c
!c       nfev is an integer output variable set to the number of
!c         calls to fcn.
!c
!c       fjac is an output n by n array which contains the
!c         orthogonal matrix q produced by the qr factorization
!c         of the final approximate jacobian.
!c
!c       ldfjac is a positive integer input variable not less than n
!c         which specifies the leading dimension of the array fjac.
!c
!c       r is an output array of length lr which contains the
!c         upper triangular matrix produced by the qr factorization
!c         of the final approximate jacobian, stored rowwise.
!c
!c       lr is a positive integer input variable not less than
!c         (n*(n+1))/2.
!c
!c       qtf is an output array of length n which contains
!c         the vector (q transpose)*fvec.
!c
!c       wa1, wa2, wa3, and wa4 are work arrays of length n.
!c
!c     subprograms called
!c
!c       user-supplied ...... fcn
!c
!c       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
!c                            qform,qrfac,r1mpyq,r1updt
!c
!c       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
!c
!c     argonne national laboratory. minpack project. march 1980.
!c     burton s. garbow, kenneth e. hillstrom, jorge j. more
!c
!c     **********
      integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,&
                      prered,p2,p5,p001,p0001,ratio,sum1,temp,xnorm,&
                      zero
      double precision dpmpar,enorm
      data one,p2,p5,p001,p0001,zero&
          /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
!c
!c     epsmch is the machine precision.
!c
      epsmch = dpmpar(1)
!c
      info = 0
      iflag = 0
      nfev = 0
!c
!c     check the input parameters for errors.
!c
      if (n .le. 0 .or. xtol .lt. zero .or. maxfev .le. 0&
         .or. ml .lt. 0 .or. mu .lt. 0 .or. factor .le. zero&
         .or. ldfjac .lt. n .or. lr .lt. (n*(n + 1))/2) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
!c
!c     evaluate the function at the starting point
!c     and calculate its norm.
!c
      iflag = 1
      call fcn(n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(n,fvec)
!c
!c     determine the number of calls to fcn needed to compute
!c     the jacobian matrix.
!c
      msum = min0(ml+mu+1,n)
!c
!c     initialize iteration counter and monitors.
!c
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
!c
!c     beginning of the outer loop.
!c
   30 continue
         jeval = .true.
!c
!c        calculate the jacobian matrix.
!c
         iflag = 2
         call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
         nfev = nfev + msum
         if (iflag .lt. 0) go to 300
!c
!c        compute the qr factorization of the jacobian.
!c
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
!c
!c        on the first iteration and if mode is 1, scale according
!c        to the norms of the columns of the initial jacobian.
!c
         if (iter .ne. 1) go to 70
         if (mode .eq. 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   40       continue
   50    continue
!c
!c        on the first iteration, calculate the norm of the scaled x
!c        and initialize the step bound delta.
!c
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   70    continue
!c
!c        form (q transpose)*fvec and store in qtf.
!c
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) .eq. zero) go to 110
            sum1 = zero
            do 90 i = j, n
               sum1 = sum1 + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum1/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
!c
!c        copy the triangular factor of the qr factorization into r.
!c
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 .lt. 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) .eq. zero) sing = .true.
  150       continue
!c
!c        accumulate the orthogonal factor in fjac.
!c
         call qform(n,n,fjac,ldfjac,wa1)
!c
!c        rescale if necessary.
!c
         if (mode .eq. 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
!c
!c        beginning of the inner loop.
!c
  180    continue
!c
!c           if requested, call fcn to enable printing of iterates.
!c
            if (nprint .le. 0) go to 190
            iflag = 0
            if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag)
            if (iflag .lt. 0) go to 300
  190       continue
!c
!c           determine the direction p.
!c
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
!c
!c           store the direction p and x + p. calculate the norm of p.
!c
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
!c
!c           on the first iteration, adjust the initial step bound.
!c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
!c
!c           evaluate the function at x + p and calculate its norm.
!c
            iflag = 1
            call fcn(n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(n,wa4)
!c
!c           compute the scaled actual reduction.
!c
            actred = -one
            if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!c
!c           compute the scaled predicted reduction.
!c
            l = 1
            do 220 i = 1, n
               sum1 = zero
               do 210 j = i, n
                  sum1 = sum1 + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum1
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
!c
!c           compute the ratio of the actual to the predicted
!c           reduction.
!c
            ratio = zero
            if (prered .gt. zero) ratio = actred/prered
!c
!c           update the step bound.
!c
            if (ratio .ge. p2) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio .ge. p5 .or. ncsuc .gt. 1) delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) .le. p2) delta = pnorm/p5
  240       continue
!c
!c           test for successful iteration.
!c
            if (ratio .lt. p0001) go to 260
!c
!c           successful iteration. update x, fvec, and their norms.
!c
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
!c
!c           determine the progress of the iteration.
!c
            nslow1 = nslow1 + 1
            if (actred .ge. p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred .ge. p2) nslow2 = 0
!c
!c           test for convergence.
!c
            if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
            if (info .ne. 0) go to 300
!c
!c           tests for termination and stringent tolerances.
!c
            if (nfev .ge. maxfev) info = 2
            if (p2*dmax1(p2*delta,pnorm) .le. epsmch*xnorm) info = 3
            if (nslow2 .eq. 10) info = 4
            if (nslow1 .eq. 20) info = 5
            if (info .ne. 0) go to 300
!c
!c           criterion for recalculating jacobian approximation
!c           by forward differences.
!c
            if (ncfail .eq. 2) go to 290
!c
!c           calculate the rank one modification to the jacobian
!c           and update qtf if necessary.
!c
            do 280 j = 1, n
               sum1 = zero
               do 270 i = 1, n
                  sum1 = sum1 + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum1 - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio .ge. p0001) qtf(j) = sum1
  280          continue
!c
!c           compute the qr factorization of the updated jacobian.
!c
            call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call r1mpyq(1,n,qtf,1,wa2,wa3)
!c
!c           end of the inner loop.
!c
            jeval = .false.
            go to 180
  290    continue
!c
!c        end of the outer loop.
!c
         go to 30
  300 continue
!c
!c     termination, either normal or user imposed.
!c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(n,x,fvec,iflag)
      return
!c
!c     last card of subroutine hybrd.
!c
end subroutine


subroutine testfun(n,x,fvec,iflag)
integer:: n,iflag,i
double precision:: x(n),fvec(n)
do i = 1,n
  fvec(i) = x(1)**3-i*x(i)+i**2
enddo
return
end


subroutine hybrd10(fcn,n,x,fvec,tol,info,wa,lwa)
  use SDRCFunctionsVar
  integer n,info,lwa
  double precision tol
  double precision x(n),fvec(n),wa(lwa)
  external fcn
!     **********
!
!    subroutine hybrd1
!
!    the purpose of hybrd1 is to find a zero of a system of
!    n nonlinear functions in n variables by a modification
!    of the powell hybrid method. this is done by using the
!    more general nonlinear equation solver hybrd. the user
!    must provide a subroutine which calculates the functions.
!    the jacobian is then calculated by a forward-difference
!    approximation.
!
!    the subroutine statement is
!
!      subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
!
!    where
!
!      fcn is the name of the user-supplied subroutine which
!        calculates the functions. fcn must be declared
!        in an external statement in the user calling
!        program, and should be written as follows.
!
!        subroutine fcn(n,x,fvec,iflag)
!        integer n,iflag
!        double precision x(n),fvec(n)
!        ----------
!        calculate the functions at x and
!        return this vector in fvec.
!        ---------
!        return
!        end
!
!        the value of iflag should not be changed by fcn unless
!        the user wants to terminate execution of hybrd1.
!        in this case set iflag to a negative integer.
!
!      n is a positive integer input variable set to the number
!        of functions and variables.
!
!      x is an array of length n. on input x must contain
!        an initial estimate of the solution vector. on output x
!        contains the final estimate of the solution vector.
!
!      fve!is an output array of length n which contains
!        the functions evaluated at the output x.
!
!      tol is a nonnegative input variable. termination occurs
!        when the algorithm estimates that the relative error
!        between x and the solution is at most tol.
!
!      info is an integer output variable. if the user has
!        terminated execution, info is set to the (negative)
!        value of iflag. see description of fcn. otherwise,
!        info is set as follows.
!
!        info = 0   improper input parameters.
!
!        info = 1   algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!        info = 2   number of calls to fcn has reached or exceeded
!                   200*(n+1).
!
!        info = 3   tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!        info = 4 or 5  iteration is not making good progress.
!
!      wa is a work array of length lwa.
!
!      lwa is a positive integer input variable not less than
!        (n*(3*n+13))/2.
!
!    subprograms called
!
!      user-supplied ...... fcn
!
!      minpack-supplied ... hybrd
!
!    argonne national laboratory. minpack project. march 1980.
!    burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!    **********
  integer index1,j,lr,maxfev,ml,mode,mu,nfev,nprint
  double precision epsfcn,factor,one,xtol,zero
  data factor,one,zero /1.0d2,1.0d0,0.0d0/
  info = 0
!
!    check the input parameters for errors.
!
  if (n .le. 0 .or. tol .lt. zero .or. lwa .lt. (n*(3*n + 13))/2) go to 20
!
!    call hybrd.
!
  maxfev = 200*(n + 1)
  xtol = tol
  ml = n - 1
  mu = n - 1
  epsfcn = zero
  mode = 2
  do 10 j = 1, n
    wa(j) = one
    10    continue
    nprint = 0
    lr = (n*(n + 1))/2
    index1 = 6*n + lr
!   print*,n,x
    call hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode,&
               factor,nprint,info,nfev,wa(index1+1),n,wa(6*n+1),lr,&
               wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
    !     if (info .eq. 5) info = 4
    20 continue
  return
end subroutine

subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
      integer n,lr
      double precision delta
      double precision r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
      integer i,j,jj,jp1,k,l
      double precision alpha,bnorm,epsmch,gnorm,one,qnorm,sgnorm,sum1,temp,zero
      double precision dpmpar,enorm
      data one,zero /1.0d0,0.0d0/
      epsmch = dpmpar(1)
      jj = (n*(n + 1))/2 + 1
      do 50 k = 1, n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum1 = zero
         if (n .lt. jp1) go to 20
         do 10 i = jp1, n
            sum1 = sum1 + r(l)*x(i)
            l = l + 1
   10       continue
   20    continue
         temp = r(jj)
         if (temp .ne. zero) go to 40
         l = j
         do 30 i = 1, j
            temp = dmax1(temp,dabs(r(l)))
            l = l + n - i
   30       continue
         temp = epsmch*temp
         if (temp .eq. zero) temp = epsmch
   40    continue
         x(j) = (qtb(j) - sum1)/temp
   50    continue
      do 60 j = 1, n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
   60    continue
      qnorm = enorm(n,wa2)
      if (qnorm .le. delta) go to 140
      l = 1
      do 80 j = 1, n
         temp = qtb(j)
         do 70 i = j, n
            wa1(i) = wa1(i) + r(l)*temp
            l = l + 1
   70       continue
         wa1(j) = wa1(j)/diag(j)
   80    continue
      gnorm = enorm(n,wa1)
      sgnorm = zero
      alpha = delta/qnorm
      if (gnorm .eq. zero) go to 120
      do 90 j = 1, n
         wa1(j) = (wa1(j)/gnorm)/diag(j)
   90    continue
      l = 1
      do 110 j = 1, n
         sum1 = zero
         do 100 i = j, n
            sum1 = sum1 + r(l)*wa1(i)
            l = l + 1
  100       continue
         wa2(j) = sum1
  110    continue
      temp = enorm(n,wa2)
      sgnorm = (gnorm/temp)/temp
      alpha = zero
      if (sgnorm .ge. delta) go to 120
      bnorm = enorm(n,qtb)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2&
            + dsqrt((temp-(delta/qnorm))**2&
                    +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp
  120 continue
      temp = (one - alpha)*dmin1(sgnorm,delta)
      do 130 j = 1, n
         x(j) = temp*wa1(j) + alpha*x(j)
  130    continue
  140 continue
      return
end subroutine

subroutine qform(m,n,q,ldq,wa)
      integer m,n,ldq
      double precision q(ldq,m),wa(m)
      integer i,j,jm1,k,l,minmn,np1
      double precision one,sum1,temp,zero
      data one,zero /1.0d0,0.0d0/
      minmn = min0(m,n)
      if (minmn .lt. 2) go to 30
      do 20 j = 2, minmn
         jm1 = j - 1
         do 10 i = 1, jm1
            q(i,j) = zero
   10       continue
   20    continue
   30 continue
      np1 = n + 1
      if (m .lt. np1) go to 60
      do 50 j = np1, m
         do 40 i = 1, m
            q(i,j) = zero
   40       continue
         q(j,j) = one
   50    continue
   60 continue
      do 120 l = 1, minmn
         k = minmn - l + 1
         do 70 i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
   70       continue
         q(k,k) = one
         if (wa(k) .eq. zero) go to 110
         do 100 j = k, m
            sum1 = zero
            do 80 i = k, m
               sum1 = sum1 + q(i,j)*wa(i)
   80          continue
            temp = sum1/wa(k)
            do 90 i = k, m
               q(i,j) = q(i,j) - temp*wa(i)
   90          continue
  100       continue
  110    continue
  120    continue
      return
end subroutine

subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      integer m,n,lda,lipvt
      integer ipvt(lipvt)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum1,temp,zero
      double precision dpmpar,enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
      epsmch = dpmpar(1)
      do 10 j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   10    continue
      minmn = min0(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
         kmax = j
         do 20 k = j, n
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k
   20       continue
         if (kmax .eq. j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm .eq. zero) go to 100
         if (a(j,j) .lt. zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
         jp1 = j + 1
         if (n .lt. jp1) go to 100
         do 90 k = jp1, n
            sum1 = zero
            do 60 i = j, m
               sum1 = sum1 + a(i,j)*a(i,k)
   60          continue
            temp = sum1/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp*temp))
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
            rdiag(k) = enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
   80       continue
   90       continue
  100    continue
         rdiag(j) = -ajnorm
  110    continue
      return
end subroutine

subroutine r1mpyq(m,n,a,lda,v,w)
      integer m,n,lda
      double precision a(lda,n),v(n),w(n)
      integer i,j,nmj,nm1
      double precision cos1,one,sin1,temp
      data one /1.0d0/
      nm1 = n - 1
      if (nm1 .lt. 1) go to 50
      do 20 nmj = 1, nm1
         j = n - nmj
         if (dabs(v(j)) .gt. one) cos1 = one/v(j)
         if (dabs(v(j)) .gt. one) sin1 = dsqrt(one-cos1*cos1)
         if (dabs(v(j)) .le. one) sin1 = v(j)
         if (dabs(v(j)) .le. one) cos1 = dsqrt(one-sin1*sin1)
         do 10 i = 1, m
            temp = cos1*a(i,j) - sin1*a(i,n)
            a(i,n) = sin1*a(i,j) + cos1*a(i,n)
            a(i,j) = temp
   10       continue
   20    continue
      do 40 j = 1, nm1
         if (dabs(w(j)) .gt. one) cos1 = one/w(j)
         if (dabs(w(j)) .gt. one) sin1 = dsqrt(one-cos1*cos1)
         if (dabs(w(j)) .le. one) sin1 = w(j)
         if (dabs(w(j)) .le. one) cos1 = dsqrt(one-sin1*sin1)
         do 30 i = 1, m
            temp = cos1*a(i,j) + sin1*a(i,n)
            a(i,n) = -sin1*a(i,j) + cos1*a(i,n)
            a(i,j) = temp
   30       continue
   40    continue
   50 continue
      return
end subroutine

subroutine r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)
      integer i,j,jj,l,nmj,nm1
      double precision cos1,cotan1,giant,one,p5,p25,sin1,tan1,tau,temp,zero
      double precision dpmpar
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
      giant = dpmpar(3)
      jj = (n*(2*m - n + 1))/2 - (m - n)
      l = jj
      do 10 i = n, m
         w(i) = s(l)
         l = l + 1
   10    continue
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 nmj = 1, nm1
         j = n - nmj
         jj = jj - (m - j + 1)
         w(j) = zero
         if (v(j) .eq. zero) go to 50
         if (dabs(v(n)) .ge. dabs(v(j))) go to 20
            cotan1 = v(n)/v(j)
            sin1 = p5/dsqrt(p25+p25*cotan1*cotan1)
            cos1 = sin1*cotan1
            tau = one
            if (dabs(cos1)*giant .gt. one) tau = one/cos1
            go to 30
   20    continue
            tan1 = v(j)/v(n)
            cos1 = p5/dsqrt(p25+p25*tan1*tan1)
            sin1 = cos1*tan1
            tau = sin1
   30    continue
         v(n) = sin1*v(j) + cos1*v(n)
         v(j) = tau
         l = jj
         do 40 i = j, m
            temp = cos1*s(l) - sin1*w(i)
            w(i) = sin1*s(l) + cos1*w(i)
            s(l) = temp
            l = l + 1
   40       continue
   50    continue
   60    continue
   70 continue
      do 80 i = 1, m
         w(i) = w(i) + v(n)*u(i)
   80    continue
      sing = .false.
      if (nm1 .lt. 1) go to 140
      do 130 j = 1, nm1
         if (w(j) .eq. zero) go to 120
         if (dabs(s(jj)) .ge. dabs(w(j))) go to 90
            cotan1 = s(jj)/w(j)
            sin1 = p5/dsqrt(p25+p25*cotan1*cotan1)
            cos1 = sin1*cotan1
            tau = one
            if (dabs(cos1)*giant .gt. one) tau = one/cos1
            go to 100
   90    continue
            tan1 = w(j)/s(jj)
            cos1 = p5/dsqrt(p25+p25*tan1*tan1)
            sin1 = cos1*tan1
            tau = sin1
  100    continue

         l = jj
         do 110 i = j, m
            temp = cos1*s(l) + sin1*w(i)
            w(i) = -sin1*s(l) + cos1*w(i)
            s(l) = temp
            l = l + 1
  110       continue

         w(j) = tau
  120    continue

         if (s(jj) .eq. zero) sing = .true.
         jj = jj + (m - j + 1)
  130    continue
  140 continue

      l = jj
      do 150 i = n, m
         s(l) = w(i)
         l = l + 1
  150    continue
      if (s(jj) .eq. zero) sing = .true.
      return

end subroutine

double precision function enorm(n,x)
      integer n
      double precision x(n)
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
            s2 = s2 + xabs*xabs
   80    continue
   90    continue
      if (s1 .eq. zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max) enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max) enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
end function

double precision function dpmpar(i)
      integer i
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      double precision dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
      data dmach(1) /2.22044604926d-16/
      data dmach(2) /2.22507385852d-308/
      data dmach(3) /1.79769313485d+308/
      dpmpar = dmach(i)
      return
end function


!! random data generating !!
subroutine normal(rnorm,mean1,sigma1)
  use SDRCFunctionsVar
  integer :: flag
  double precision :: u1, u2, y1, y2, rnorm, mean1, sigma1
  save flag
  data flag /0/
  call RANDOM_NUMBER(u1)
  call RANDOM_NUMBER(u2)
  if (flag.eq.0) then
    y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
    rnorm = mean1 + sigma1*y1
    flag = 1
  else
    y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
    rnorm = mean1 + sigma1*y2
    flag = 0
  end if
end subroutine normal

subroutine init_random_seed()
  integer :: i, m, clock
  integer, dimension(:), allocatable :: seed

  call random_seed(size = m)
  allocate(seed(m))
  call system_clock(count=clock)
  seed = clock + 37 * (/ (i - 1, i = 1, m) /)
  call random_seed(put = seed)
  deallocate(seed)
end subroutine init_random_seed


subroutine random_exp(fn_val,para,m)
!  The function random_exp() returns a exponential distributed pseudo-random
!  number with mean para, that is the pdf given by \lambda*exp(-\lambda*x)
!  The sample size is m
  use SDRCFunctionsVar
  double precision, dimension(m), INTENT(OUT) :: fn_val
  double precision :: u
  INTEGER :: i
  !     Local variables
  !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
  call init_random_seed()
  DO i = 1,m
    call RANDOM_NUMBER(u)
    fn_val(i) = -log(1-u)/para
  ENDDO
END subroutine

 subroutine matInv(am,cm,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
  implicit none 
  integer :: n
  double precision,dimension(n,n):: am
  double precision,dimension(n,n), intent(out) :: cm
  double precision,dimension(n,n):: Lm, Um
  double precision,dimension(n):: bm, dm, xm
  double precision :: coeff
  integer :: i, j, k

  ! step 0: initialization for matrices L and U and b
  ! Fortran 90/95 allows such operations on matrices
  Lm=0.0
  Um=0.0
  bm=0.0

  ! step 1: forward elimination
  do k=1, n-1
    do i=k+1,n
      coeff=am(i,k)/am(k,k)
      Lm(i,k) = coeff
      do j=k+1,n
        am(i,j) = am(i,j)-coeff*am(k,j)
      enddo
    enddo
  enddo

  ! Step 2: prepare L and U matrices 
  ! L matrix is a matrix of the elimination coefficient
  ! + the diagonal elements are 1.0
  do i=1,n
    Lm(i,i) = 1.0
  enddo
  ! U matrix is the upper triangular part of A
  do j=1,n
    do i=1,j
      Um(i,j) = am(i,j)
    enddo
  enddo

  ! Step 3: compute columns of the inverse matrix C
  do k=1,n
    bm(k)=1.0
    dm(1) = bm(1)
    ! Step 3a: Solve Ld=b using the forward substitution
    do i=2,n
      dm(i)=bm(i)
      do j=1,i-1
        dm(i) = dm(i) - Lm(i,j)*dm(j)
      enddo
    enddo
    ! Step 3b: Solve Ux=d using the back substitution
    xm(n)=dm(n)/Um(n,n)
    do i = n-1,1,-1
      xm(i) = dm(i)
      do j=n,i+1,-1
        xm(i)=xm(i)-Um(i,j)*xm(j)
      enddo
      xm(i) = xm(i)/um(i,i)
    enddo
    ! Step 3c: fill the solutions x(n) into column k of C
    do i=1,n
      cm(i,k) = xm(i)
    enddo
    bm(k)=0.0
  enddo
end subroutine matInv


subroutine kernel(xVal,m,kernelVal1,typeKert,sized)
! this is the kernel function. m is the length of vector xVal, d is the dimension of each obsevation vector xVal
! typeKer is type of kernels, including: uniform, quartic, derivative of quartic, normal
  use SDRCFunctionsVar
  integer :: i,j,k, t,m,sized
  double precision, dimension(m,sized) :: xVal
  double precision :: valTemp, indicVal, c, valTemp1
  double precision, dimension(m), intent(out):: kernelVal1
  character (len=4) :: typeKert
  c = 0.9375
  select case (typeKert)
  !case ('uniform') ! this is for uniform kernel.
  ! where (abs(xVal) <= 1)
  !   kernelVal = 0.5
  ! else where
  !   kernelVal = 0.0
  ! end where
  case ('quar') ! this is quaratic kernel
    do i = 1,m
      kernelVal1(i) = 1.0
      do j = 1,sized
        if (xVal(i,j)<= 1.0 .AND. xVal(i,j)>= -1.0) then
          valTemp =  0.9375*(1.0-xVal(i,j)*xVal(i,j))*(1.0-xVal(i,j)*xVal(i,j))
        else
          valTemp = 0.0
        end if
        kernelVal1(i) = kernelVal1(i)*valTemp
      enddo
    enddo
  case ('norm') ! this use pdf of normal distribution as kernel
    do i = 1,m
      kernelVal1(i) = 1.0
      do j = 1,sized
        valTemp = 1/sqrt(2.0*pi)*exp(-xVal(i,j)**2/2.0)
        kernelVal1(i) = kernelVal1(i)*valTemp
      enddo
    enddo
  !case ('Ep') ! this use Epanechnikov kernel
  ! where (abs(xVal) <= 1)
  !   kernelVal = (1-xVal**2)*3/4
  ! else where
  !   kernelVal = 0.0
  ! end where
  end select
end subroutine

subroutine kernel2(xVal,m,kernelVal2)
! this is the kernel function. n is the length of vector xVal, d is the dimension of each obsevation vector xVal
! typeKer is type of kernels, including: uniform, quartic, derivative of quartic, normal
  use SDRCFunctionsVar
  integer :: i,j, k, t, m
  double precision, dimension(m,d) :: xVal
  double precision :: valTemp, indicVal, c, valTemp1
  double precision, dimension(m,d), intent(out):: kernelVal2
  c = 0.9375
! select case (typeKer1)
  !case ('uniform') ! this is for uniform kernel.
  ! where (abs(xVal) <= 1)
  !   kernelVal = 0.5
  ! else where
  !   kernelVal = 0.0
  ! end where
! case ('quar1') ! this is first derivative of quaratic kernel
    kernelVal2 = 1.0
    do i = 1,m
      do j = 1,d
        do t = 1,d
          if (xVal(i,j)<= 1.0 .AND. xVal(i,j)>= -1.0) then
            if (t == j) then
              valTemp = 4.0*c*xVal(i,t)*(xVal(i,t)*xVal(i,t)-1.0)
            else
              valTemp = c*(1.0-xVal(i,t)*xVal(i,t))*(1.0-xVal(i,t)*xVal(i,t))
            end if
          else
            valTemp = 0.0
          end if
          kernelVal2(i,j) = kernelVal2(i,j)*valTemp
        enddo
      enddo
    enddo
  !case ('normal') ! this use pdf of normal distribution as kernel
  ! kernelVal = exp(-xVal**2/2)/sqrt(2*pi)
  !case ('Ep') ! this use Epanechnikov kernel
  ! where (abs(xVal) <= 1)
  !   kernelVal = (1-xVal**2)*3/4
  ! else where
  !   kernelVal = 0.0
  ! end where
! end select
end subroutine

!!!! functions for algorithm
! score function


