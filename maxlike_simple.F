c $Log: maxlike_simple.F,v $
c Revision 1.1.1.1  2009/03/09 12:48:28  rbb24
c
c new import of diffusion code
c
c
c Revision 1.4  2008/11/18 14:58:47  hummer
c added ranlux as random number generator
c
c added option SPREAD to allow a Gaussian spread of the initial condition
c (which is otherwise assumed to be a "delta" function); var0 is
c optimized;
c
c added option LAG to vary the zero time after which diffusion
c begins (normally _TIME0_=0; values larger than zero would be expected
c in cases where inertia is important; values smaller than zero would
c be expected for cage diffusion);
c
c changed move scheme for equilibrium probabilities and rates to being
c (more strictly) linear on a log scale
c
c Revision 1.3  2007/11/30 21:55:20  hummer
c *** empty log message ***
c
c Revision 1.2  2007/10/24 13:25:20  hummer
c now catches possibility that pt(j)==0 && mm==0.
c
c Revision 1.4  2004/09/21 17:33:27  hummer
c *** empty log message ***
c
c Revision 1.2  2004/03/01 18:32:43  hummer
c include periodic linear chain
c
c Revision 1.1  2004/02/17 15:37:15  hummer
c Initial revision
c
c $Id: maxlike_simple.F,v 1.1.1.1 2009/03/09 12:48:28 rbb24 Exp $

#ifndef _TIME_                        /* time interval */
#  define _TIME_ 20.0
#endif
#ifndef _N_                           /* dimension */
#  define _N_ 32
#endif
#ifndef _TEMP0_                       /* initial temperature */
#  define _TEMP0_ 1.0
#endif
#ifndef _TEMP1_                       /* final temperature */
#  define _TEMP1_ 0.0
#endif
#ifndef _NTIM_                        /* number of transition matrices */
#  define _NTIM_ 1
#endif
#ifndef _DR_                         /* MC jump width rate coefficient */
#  define _DR_ 0.01
#endif
#ifndef _DP_                         /* MC jump width equilibrium probability*/
#  define _DP_ 0.01
#endif
#ifndef _NMC_                        /* number of Monte Carlo steps */
#  define _NMC_ 1000
#endif
#ifndef _DIFF_
#  define _DIFF_ 0.02                 /* diffusion coefficient */
#endif
#ifndef _VAR0_                        /* initial spread var0 */
#  define _VAR0_ 1.d-1  
#endif
#ifndef _TIME0_                       /* initial lag shift time0 */
#  define _TIME0_ 0.0d0
#endif
#ifndef _DVAR0_                       /* MC move width (var0) */
#  define _DVAR0_ 1.d0
#endif
#ifndef _DTIME0_                      /* MC move width (time0) */
#  define _DTIME0_ 0.1d0
#endif

      program maximum_likelihood_MC_diffusion
c **************************************************************************
c * fit rate matrix to observed transition matrix uxing maximum likelihood *
c * uses simulated annealing                                               *
c * Gerhard Hummer, February 2004, National Institutes of Health           *
c **************************************************************************

      implicit real*8 ( a-h, o-z )
      parameter ( TEMP0 = _TEMP0_ , TEMP1 = _TEMP1_ )
      parameter ( n = _N_ )
      parameter ( NMC = _NMC_ )      ! number of MC steps
      real*8 r ( n * ( n - 1 ) / 2 ) ! free coefficients of rate matrix
      real*8 peq ( n )               ! equilibrium probability
      real*8 rsave ( n * ( n - 1 ) / 2 ) ! saved rate matrix in Monte Carlo
      real*8 peqsave ( n )           ! saved equilibrium distribution in MC
      real*8 rmax ( n * ( n - 1 ) / 2 ) ! max-like rate matrix in Monte Carlo
      real*8 peqmax ( n )           ! max-like equilibrium distribution in MC
      real*8 max_log_like            ! maximum log-likelihood
      real*8 rmat ( n, n )           ! rate matrix
      real*8 rsym ( n, n )           ! symmetrized rate matrix
      real*8 d ( n )                 ! eigenvectors
      real*8 v ( n, n )              ! eigenvalues
      real*8 fv1 ( n ), fv2 ( n )    ! temporary
      integer ind ( n )              ! index vector
      parameter ( ntim = _NTIM_ )    ! number of transition matrices read
      real*8 mm ( n, n, ntim )       ! transition matrices
      real*8 m                       ! normalization constant
      real*8 log_like                ! log likelihood
      real*8 log_like_try            ! trial log likelihood
      parameter ( dr = _DR_ )        ! MC jump width rate coefficient
      parameter ( dp = _DP_ )        ! MC jump width equilibrium probability
c
      integer iaccr                  ! number of accepted MC moves (rate)
      integer iaccp                  ! number of accepted MC moves (probabil)
      integer iaccsp                 ! number of accepted MC moves (spread)
      logical accept                 ! logical (accept)
      character*255 datfile          ! file with transition matrices
c
      real*8 piconst                 ! Pi
      parameter ( piconst = 3.14159265358979d0 )
      real*8 diff ( n )              ! diffusion coefficient (assum. 2*Pi dist)
      real*8 var0                    ! variance of t=0 spread in diffusion
                                     ! used in SPREAD is defined
      real*8 time0                   ! shift in lag time used if LAG is defined
      integer iran
c
      write ( 6, 1 )
 1    format (
     $ '************************************************************'/,
     $ 'Fit rate matrix to observed transition matrix'/,
     $ 'uxing maximum likelihood;  Uses simulated annealing'/,
     $ 'Gerhard Hummer, February 2004, NIH'/,
     $ 'NOTE: input transition matrix has transitions i->j in ROW i'/,
     $ '************************************************************' )
c
      call getarg ( 1, datfile )
      write ( 6, 2 ) datfile ( 1 : lastc ( datfile ) )
 2    format ( 'file with transition matrices:', a )
c output info
      write ( 6, 3 ) _N_, _TIME_, _TEMP0_, _TEMP1_, _NTIM_, _DR_,
     $ _DP_, _NMC_
 3    format (
     $ 'dimension =                     ', i14/,
     $ 'time step between transitions = ', e14.6/,
     $ 'initial temperature =           ', e14.6/,
     $ 'final temperature =             ', e14.6/,
     $ 'number of transition matrices = ', i14/,
     $ 'MC move width (rate) =          ', e14.6/,
     $ 'MC move width (equ prob) =      ', e14.6/,
     $ 'number of MC steps =            ', i14 )
      time0 = _TIME0_
      var0 = _VAR0_
      dvar0 = _DVAR0_
      dtime0 = _DTIME0_
c open output file
      open ( 66, file = 'maximum_likelihood_MC_diffusion.log',
     $ form = 'formatted', status = 'unknown' )
c initiate new configuration
        ppp = 0.d0
        do i = 1, n
          peq ( i ) = ranx ( iran ) / dble ( n ) + 0.001
          ppp = ppp + peq ( i )
        enddo
        do i = 1, n
          peq ( i ) = peq ( i ) / ppp
        enddo
        do i = 1, n * ( n - 1 ) / 2
          r ( i ) = .05 * ranx ( iran ) + 0.0001
        enddo
c set certain rate coefficients to zero (here: linear chain with
c                                  periodic boundary conditions)
#ifdef LINEAR
      l = 0
      ll = 0
      do i = 2, n
        do j = 1, i - 1
          l = l + 1
          ij = i - j
          if ( ij .gt. 1 ) then
            r ( l ) = 0.
          else
            ll = ll + 1
            write ( 0, * ) i, j
          endif
        enddo
      enddo
#endif
      ll = 0
      do l = 1, n * ( n - 1 ) / 2
        if ( r ( l ) .gt. 0. ) then
          ll = ll + 1
        endif
      enddo
      write ( 6, 4 ) ll
 4    format ( 'number of variable coefficients in rate matrix =', i8 )
      write ( 6, 5 ) n - 1
 5    format ( 'number of variable coefficients in p-eq =', i8 )

c read in transition matrix
      open ( 10, file = datfile,
     $ form = 'formatted', status = 'old' )
c
      do itim = 1, ntim
        do i = 1, n
          read ( 10, * ) ( mm ( j, i, itim ), j = 1, n )
        enddo
        do i = 1, n
          do j = 1, n
            write ( 99, '(100i12)' )
     $       i, j, itim, int ( mm ( j, i, itim ) ), 
     $       int ( mm ( i, j, itim ) )
          enddo
        enddo
      enddo
c
      close ( 10 )
c
#ifdef LINEAR
#ifdef _DELTA_
      delta = _DELTA_
      write ( 6, '(''grid mesh width = '', 1pe14.6 )' ) delta
      if ( iargc ( ) .le. 1 ) then
c set up equilibrium probability by applying detailed balance to tmat
        do j = 1, n
          m = 0.d0
          do i = 1, n
            m = m + mm ( i, j, ntim )
          enddo
          if ( m .gt. 0. ) then
            do i = 1, n
              rmat ( i, j ) = mm ( i, j, ntim ) / m
            enddo
          else
            do i = 1, n
              rmat ( i, j ) = 0.0
            enddo
          endif
        enddo
        do j = 1, n
          do i = 1, n
            if ( rmat ( j, i ) .gt. 0. ) then
              rsym ( j, i ) = rmat ( j, i )
     $         * sqrt ( rmat ( i, j ) / rmat ( j, i ) )
            else
              rsym ( j, i ) = rmat ( j, i )
            endif
          enddo
        enddo
        call rs ( n, n, rsym, d, 1, v, fv1, fv2, ierr )
        call indexx ( n, d, ind )
        ppp = 0.d0
        do i = 1, n
          ppp = ppp + max ( abs ( v ( i, ind ( n ) ) ), 1.d-6 )
        enddo
        do i = 1, n
          peq ( i ) = max ( abs ( v ( i, ind ( n ) ) ), 1.d-6 ) / ppp
        enddo
c set up rates using initial D0 and equilibrium probabilities
c
        l = 0
        do i = 2, n
          do j = 1, i - 1
            l = l + 1
            ij = i - j
            if ( ij .gt. 1 ) then
              r ( l ) = 0.
            else
              r ( l ) =  ( _DIFF_ / delta ** 2 ) *
     $         sqrt ( peq ( i ) / peq ( j ) )
     $         * ( 1. + ( ranx ( iran ) - 0.5 ) / 10. )
            endif
          enddo
        enddo
      endif
c construct rate matrix...
      l = 0
      do i = 2, n
        do j = 1, i - 1
          l = l + 1
          rmat ( i, j ) = r ( l )
          rmat ( j, i ) = r ( l ) * peq ( j ) / peq ( i )
        enddo
      enddo
c calculate intitial local diffusion coefficients
      do i = 1, n
        if ( i .lt. n ) then
          diff ( i ) = rmat ( i + 1, i )
     $     * sqrt ( peq ( i ) / peq ( i + 1 ) ) *
     $     delta ** 2
        else
          diff ( i ) = 0.
        endif
      enddo
      do i = 1, n
        write ( 6, * ) i, peq ( i ), diff ( i )
      enddo
#endif

C calculate initial log likelihood
      call log_likelihood ( r, peq, mm, log_like, s, delta, var0,
     $ time0 )
c
      write ( 6, * ) 'initial log-like = ', log_like
c maximum likelihood
      max_log_like = log_like   ! initialize maximum
      do i = 1, n * ( n - 1 ) / 2
        rmax ( i ) = r ( i )
      enddo
      do i = 1, n
        peqmax ( i ) = peq ( i )
      enddo
c
      iran = 123456789          ! random number seed
      iaccr = 0                 ! accepted R moves
      iaccp = 0                 ! accepted peq moves
      iaccsp = 0                ! accepted spread moves
      iacclag = 0               ! accepted lag time moves
C MAIN MC LOOP
C each loop does one R move and one Peq move
      do imc = 0, NMC
c linearly interpolate temperature for simulated annealing
        temp = ( ( NMC - imc ) * TEMP0 + imc * TEMP1 ) / NMC
c move r
 111    j = ranx ( iran ) * ( n * ( n - 1 ) ) / 2 + 1
        j = min ( j, ( n * ( n - 1 ) ) / 2 )
        if ( r ( j ) .le. 0.d0 ) goto 111 ! only modify rate coefficients > 0
        rsave ( j ) = r ( j )
cGH        ddr = dr * ( ranx ( iran ) - 0.5 )
cGH        if ( r ( j ) + ddr .gt. 0. ) then
cGH          r ( j ) = r ( j ) + ddr
cGH          r ( j ) = r ( j ) * abs ( 1.d0 + ddr )
c move r(j) on a log scale
        ddr = dr * ( ranx ( iran ) - 0.5 )
        !print *,ddr
        r ( j ) = r ( j ) * exp ( ddr )
          call log_likelihood ( r, peq, mm, log_like_try, s, delta,
     $     var0, time0 )
c Metropolis
          accept = log_like_try .gt. log_like
          if ( .not. accept ) then
            ppp = exp ( ( log_like_try - log_like ) / temp )
            rrr = ranx ( iran )
            accept = ppp .gt. rrr
          endif
          if ( accept ) then
            log_like = log_like_try
            iaccr = iaccr + 1
c maximum likelihood
            if ( max_log_like .lt. log_like ) then ! update maximum
              max_log_like = log_like
              do i = 1, n * ( n - 1 ) / 2
                rmax ( i ) = r ( i )
              enddo
              do i = 1, n
                peqmax ( i ) = peq ( i )
              enddo
            endif
          else
            r ( j ) = rsave ( j )
          endif
cGH        endif
c move p
        do i = 1, n
          peqsave ( i ) = peq ( i )
        enddo
        j = ranx ( iran ) * n + 1
        j = min ( j, n )
cGH        ddp = log ( dp ) * ( ranx ( iran ) - 0.5 )
c move p(j) on a log scale
        ddp = dp * ( ranx ( iran ) - 0.5 )
        peq ( j ) = peq ( j ) * exp ( ddp )
        ppp = 0.
        do i = 1, n
          ppp = ppp + peq ( i )
        enddo
        do i = 1, n
          peq ( i ) = peq ( i ) / ppp
        enddo
        call log_likelihood ( r, peq, mm, log_like_try, s, delta,
     $   var0, time0 )
c Metropolis
        accept = log_like_try .gt. log_like
        if ( .not. accept ) then
          ppp = exp ( ( log_like_try - log_like ) / temp )
          rrr = ranx ( iran )
          accept = ppp .gt. rrr
        endif
        if ( accept ) then
          log_like = log_like_try
          iaccp = iaccp + 1
c maximum likelihood
          if ( max_log_like .lt. log_like ) then ! update maximum
            max_log_like = log_like
            do i = 1, n * ( n - 1 ) / 2
              rmax ( i ) = r ( i )
            enddo
            do i = 1, n
              peqmax ( i ) = peq ( i )
            enddo
          endif
        else
          do i = 1, n
            peq ( i ) = peqsave ( i )
          enddo
        endif
c print out P and D to log file ...
        if ( mod ( imc, 100 ) .eq. 0 ) then
c diagonalize
          l = 0
          do i = 2, n
            do j = 1, i - 1
              l = l + 1
              rmat ( i, j ) = r ( l )
              rmat ( j, i ) = r ( l ) * peq ( j ) / peq ( i )
            enddo
          enddo
          do i = 1, n
            rmat ( i, i ) = 0.
            do j = 1, n
              if ( i .ne. j ) then
                rmat ( i, i ) = rmat ( i, i ) - rmat ( j, i )
              endif
            enddo
          enddo
          do i = 1, n
            do j = 1, n
              if ( rmat ( j, i ) .gt. 0. ) then
                rsym ( j, i ) = rmat ( j, i )
     $           * sqrt ( rmat ( i, j ) / rmat ( j, i ) )
              else
                rsym ( j, i ) = rmat ( j, i )
              endif
            enddo
          enddo
          call rs ( n, n, rsym, d, 1, v, fv1, fv2, ierr )
          call indexx ( n, d, ind )
          if ( mod ( imc, 1000 ) .eq. 0 ) then          
            write ( 0,  '(i8,x,e15.8,x,100(e12.6,x))' ) imc, -log_like,
     $       -1./d(ind(n-1)),
     $       r ( 1 ), ( peq ( i ), i = 1, min ( 5, n ) )
            write ( 6,  '(i8,x,e15.8,x,100(e12.6,x))' ) imc, -log_like,
     $       -1./d(ind(n-1)),
     $       r ( 1 ), ( peq ( i ), i = 1, min ( 5, n ) ),
     $       var0
          endif
#ifdef LINEAR
          do i = 1, n - 1
            diff ( i ) = rmat ( i + 1, i )
     $       * sqrt ( peq ( i ) / peq ( i + 1 ) ) *
     $       delta ** 2
          enddo
          diff ( n ) = 0.
#endif
c
          write ( 66, '(i8,1000(1pe14.6))' ) imc, -log_like,
     $     -1./d(ind(n-1)),
     $     ( peq ( i ), i = 1, n ),
     $     ( diff ( i ), i = 1, n )
#else
          write ( 66,  '(i8,100(1pe14.6))' ) imc, -log_like,
     $     -1./d(ind(n-1)),
     $     r ( 1 ), ( peq ( i ), i = 1, min ( 10, n ) )
#endif
        endif
      enddo
      ! END OF MAIN MC LOOP
      write ( 6, * ) 'number of accepted rate moves = ',iaccr
      write ( 6, * ) 'number of accepted peq moves =  ',iaccp
      write ( 6, * ) 'maximum log-likelihood = ', max_log_like
      write ( 6, * ) 'final log-likelihood = ', log_like
      write ( 6, * ) 'expected log-likelihood = ', s
c
      write ( 6, '(10(1pe16.9))' )( peqmax ( i ), i = 1, n )
      write ( 6, '(10(1pe16.9))' )( rmax ( i ),
     $ i = 1, n * ( n - 1 ) / 2 )
c
      write ( 6, '(10(1pe16.9))' )( peq ( i ), i = 1, n )
      write ( 6, '(10(1pe16.9))' )( r ( i ), i = 1, n * ( n - 1 ) / 2 )
      write ( 6, '(10(1pe16.9))' )( r ( i ) * peq ( i ) /
     $ peq ( i + 1 ),
     $ i = 1, n - 1 )
c diagonalize
      l = 0
      do i = 2, n
        do j = 1, i - 1
          l = l + 1
          rmat ( i, j ) = r ( l )
          rmat ( j, i ) = r ( l ) * peq ( j ) / peq ( i )
        enddo
      enddo
      do i = 1, n
        rmat ( i, i ) = 0.
        do j = 1, n
          if ( i .ne. j ) then
            rmat ( i, i ) = rmat ( i, i ) - rmat ( j, i )
          endif
        enddo
      enddo
      do i = 1, n
        do j = 1, n
          if ( rmat ( j, i ) .gt. 0. ) then
            rsym ( j, i ) = rmat ( j, i )
     $       * sqrt ( rmat ( i, j ) / rmat ( j, i ) )
          else
            rsym ( j, i ) = rmat ( j, i )
          endif
        enddo
      enddo
      call rs ( n, n, rsym, d, 1, v, fv1, fv2, ierr )
      call indexx ( n, d, ind )
c save rate matrix
      open ( 10, file = 'rate.dat', form = 'formatted',
     $ status = 'unknown' )
      write ( 10, '(i6)' ) n
      write ( 10, '(5(1pe15.7))' )
     $ ( ( rmat ( i, j ), j = 1, n ), i = 1, n )
      close ( 10 )
c save eigensystem
      open ( 10, file = 'rate_eigen.dat', form = 'formatted',
     $ status = 'unknown' )
      write ( 10, '(i6)' ) n
      write ( 10, '(5(1pe15.7))' ) ( d ( ind ( i ) ), i = n, 1, -1 )
      write ( 10, '(5(1pe15.7))' )
     $ ( ( v ( i, ind ( j ) ), j = n, 1, -1 ), i = 1, n )
      close ( 10 )
c save rate + peq
      open ( 10, file = 'maximum_likelihood_MC_diffusion.save',
     $ form = 'formatted', status = 'unknown' )
      write ( 10, '(4(1pe20.12))' ) ( peq ( i ), i = 1, n  )
      write ( 10, '(4(1pe20.12))' )
     $ ( r ( i ), i = 1, n * ( n - 1 ) / 2 )
#ifdef SPREAD
      write ( 10, '(1pe20.12)' ) var0
#endif
#ifdef LAG
      write ( 10, '(1pe20.12)' ) time0
#endif
      close ( 10 )
c
      close ( 66 )
c
      end

      subroutine log_likelihood ( r, peq, mm, log_like, s, delta,
     $ var0, time0 )
      implicit real*8 ( a-h, o-z )
      parameter ( n = _N_ )
      real*8 r ( n * ( n - 1 ) / 2 ) ! free coefficients of rate matrix
      real*8 peq ( n )               ! equilibrium probability
      real*8 p ( n )                 ! scaled probabilities at time 0
      parameter ( ntim = _NTIM_ )    ! number of transition matrices read
      real*8 mm ( n, n, ntim )       ! transition matrices
      real*8 m                       ! normalization constant
      real*8 rmat ( n, n )           ! rate matrix
      real*8 rsym ( n, n )           ! symmetrized rate matrix
      real*8 d ( n )                 ! eigenvectors
      real*8 v ( n, n )              ! eigenvalues
      real*8 pt ( n )                ! probability at time t
      real*8 fv1 ( n ), fv2 ( n )    ! temporary
      real*8 pdotv ( n )             ! dot products with eigenvectors
c
      real*8 log_like                ! log likelihood
#ifdef DIFFSTIFF
      real*8 diff ( n )              ! diffusion coefficient
#endif
#ifdef SPREAD
      real*8 pinit ( n )
#endif
      real*8 var0
      real*8 time0
c
c normalize equilibrium probabilities
      ppp = 0.
      do i = 1, n
        ppp = ppp + peq ( i )
      enddo
      do i = 1, n
        peq ( i ) = peq ( i ) / ppp
      enddo
c construct rate matrix
      l = 0
      do i = 2, n
        do j = 1, i - 1
          l = l + 1
          rmat ( i, j ) = r ( l )
          rmat ( j, i ) = r ( l ) * peq ( j ) / peq ( i )
        enddo
      enddo
      do i = 1, n
        rmat ( i, i ) = 0.
        do j = 1, n
          if ( i .ne. j ) then
            rmat ( i, i ) = rmat ( i, i ) - rmat ( j, i )
          endif
        enddo
      enddo
c
c diagonalize symmetrized matrix
      do i = 1, n
        do j = 1, n
          if ( rmat ( j, i ) .gt. 0. ) then
            rsym ( j, i ) = rmat ( j, i )
     $       * sqrt ( rmat ( i, j ) / rmat ( j, i ) )
          else
            rsym ( j, i ) = rmat ( j, i )
          endif
        enddo
      enddo
c
c      write ( 0, * ) 'before rs'
c
      call rs ( n, n, rsym, d, 1, v, fv1, fv2, ierr )
c        do i = 1, n
c          write ( 0, * ) i, d ( i )
c        enddo
c
c      write ( 0, * ) 'after rs'
c
c      call jacobi ( rsym, n, n, d, v, nrot )

      log_like = 0.
      s = 0.                              ! expected likelihood
      do i = 1, n
#ifdef SPREAD /* spread data according to Gaussian */
        psum = 0.d0
        do j = 1, n
          p ( j ) = exp ( - 0.5 * ( j - i ) ** 2 * delta ** 2 / var0 )
     $     * peq ( j )          ! Gaussian weighted peq
          psum = psum + p ( j )
        enddo
        do j = 1, n
          p ( j ) = p ( j ) / ( psum * sqrt ( peq ( j ) ) )
c          write ( 0, * ) j * delta, p ( j ) * sqrt ( peq ( j ) )
        enddo
c
        do j = 1, n                       ! dot product with eigenvectors
          do jj = 1, n
            pdotv ( j ) = p ( jj ) * v ( jj, j )
          enddo
        enddo
#else
        do j = 1, n
          p ( j ) = 0.
        enddo
        p ( i ) = 1. / sqrt ( peq ( i ) ) ! initial vector (time 0)
        do j = 1, n                       ! dot product with eigenvectors
          pdotv ( j ) = p ( i ) * v ( i, j )
        enddo
#endif
        do itim = 1, ntim
          time = itim * _TIME_ - time0
c          write ( 0, * ) time
c probabilities at time t
c pt(j) = sqrt ( peq ( j ) ) * sum_l exp(d(l)*t) * pdotv(l) * v(j,l)
          do j = 1, n
            pt ( j ) = 0.
            do l = 1, n
              pt ( j ) = pt ( j )
     $         + exp ( d ( l ) * time ) * pdotv ( l ) * v ( j, l )
            enddo
ccc            pt ( j ) = pt ( j ) * sqrt ( peq ( j ) )
ccc !!! absolute values !!!
c            write ( 0, * ) j, pt ( j )
            pt ( j ) = abs ( pt ( j ) * sqrt ( peq ( j ) ) )
ccc
          enddo
c log likelihood: number of observations * log ( probability )
          do j = 1, n
            if ( mm ( j, i, itim ) .gt. 0 ) then
              log_like = log_like +
     $         mm ( j, i, itim ) * log ( pt ( j ) )
            endif
          enddo
c expected likelihood
          m = 0.d0
          do j = 1, n
            m = m + mm ( j, i, itim )
          enddo
          do j = 1, n
            s = s + m * pt ( j ) * log ( pt ( j ) )
          enddo
        enddo
c
      enddo
c
#ifdef DIFFSTIFF
#ifdef PERIODIC
      do i = 1, n
#else
      do i = 1, n - 1
#endif
        if ( i .lt. n ) then
          diff ( i ) = rmat ( i + 1, i )
     $     * sqrt ( peq ( i ) / peq ( i + 1 ) ) *
     $     delta ** 2
        else
          diff ( i ) = rmat ( 1, i )
     $     * sqrt ( peq ( i ) / peq ( 1 ) ) *
     $     delta ** 2
        endif
      enddo
      xlogpost = 0.d0
#ifdef PERIODIC
      do i = 1, n
#else
      do i = 1, n - 2
#endif
        if ( i .lt. n ) then
          xlogpost = xlogpost
     $     - 0.5d0 * ( diff ( i ) - diff ( i + 1 ) ) ** 2
     $     / DIFFSTIFF ** 2
#ifdef PERIODIC
        else
          xlogpost = xlogpost
     $     - 0.5d0 * ( diff ( i ) - diff ( 1 ) ) ** 2
     $     / DIFFSTIFF ** 2
#endif
        endif
      enddo
      log_like = log_like + xlogpost
#endif /* DIFFSTIFF */
c
      return
      end

      SUBROUTINE indexx(n,arr,indx)
      implicit real*8 ( a-h, o-z )
      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find the eigenvalues and eigenvectors (if desired)
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        matz  is an integer variable set equal to zero if
c        only eigenvalues are desired.  otherwise it is set to
c        any non-zero integer for both eigenvalues and eigenvectors.
c
c     on output
c
c        w  contains the eigenvalues in ascending order.
c
c        z  contains the eigenvectors if matz is not zero.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat
c           and tql2.  the normal completion code is zero.
c
c        fv1  and  fv2  are temporary storage arrays.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      if (n .le. nm) go to 10
      ierr = 10 * n
      go to 50
c
   10 if (matz .ne. 0) go to 20
c     .......... find eigenvalues only ..........
      call  tred1(nm,n,a,w,fv1,fv2)
*  tqlrat encounters catastrophic underflow on the Vax
*     call  tqlrat(n,w,fv2,ierr)
      call  tql1(n,w,fv1,ierr)
      go to 50
c     .......... find both eigenvalues and eigenvectors ..........
   20 call  tred2(nm,n,a,w,fv1,z)
      call  tql2(nm,n,w,fv1,z,ierr)
   50 return
      end
      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
      subroutine tred2(nm,n,a,d,e,z)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
      subroutine tql1(n,d,e,ierr)
c
      integer i,j,l,m,n,ii,l1,l2,mml,ierr
      double precision d(n),e(n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql1,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the ql method.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        e has been destroyed.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 290 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 210
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  210    p = d(l) + f
c     .......... order eigenvalues ..........
         if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
         do 230 ii = 2, l
            i = l + 2 - ii
            if (p .ge. d(i-1)) go to 270
            d(i) = d(i-1)
  230    continue
c
  250    i = 1
  270    d(i) = p
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end

      function lastc(line)
c ***********************************************************************
c * Function returning the index of the last non-blank character in a   *
c * string                                                              *
c *                                                                     *
c * Input:                                                              *
c *                                                                     *
c *   line:   character string                                          *
c *                                                                     *
c * Output:                                                             *
c *                                                                     *
c *   lastc:  index of the last non-blank character in line             *
c *                                                                     *
c * Author:                Gerhard Hummer, Los Alamos National Lab      *
c * Date:                  09-OCT-98                                    *
c ***********************************************************************

      implicit none
      character*(*) line
      integer lastc             ! index of last character in a string
      lastc=len(line)
 1    if(line(lastc:lastc).ne.' ')goto 2
      lastc=lastc-1
      if(lastc.ge.1)goto 1
 2    lastc=max(lastc,1)
      return
      end

c *** luxury random number generator ***

      function ranx(iran)
      implicit real*8 (a-h,o-z)
      parameter(f=1.d0/2.d0**24)
c      iran=iran*69069+1
      iran=iran*1664525+1013904223
      ranx=ishft(iran,-8)*f
      return
      end

C      function ranx ( iran )
C      implicit none
C      integer n, isw, iran
C      real ranx
C      parameter ( n = 2**12 )
C      real x ( n )
C      data isw / 0 /
C      save isw, x
C      if ( isw .gt. 0 ) then
C        if ( isw .le. n ) then
C          ranx = x ( isw )
C          isw  = isw + 1
C        else
C          call RANLUX ( x, n )
C          ranx = x ( 1 )
C          isw  = 2
C        endif
C      else
Cc        call rluxgo ( 3, iran, 0, 0 )
C        call rluxgo ( 3, iran, 0, 0 )
C        call RANLUX ( x, n )
C        ranx = x ( 1 )
C        isw  = 2
C      endif
C      return
C      end

      SUBROUTINE RANLUX(RVEC,LENV)                                          0084
c                                                                           0085
C     Subtract-and-borrow random number generator proposed by               0086
C     Marsaglia and Zaman, implemented by F. James with the name            0087
C     RCARRY in 1991, and later improved by Martin Luescher                 0088
C     in 1993 to produce "Luxury Pseudorandom Numbers".                     0089
c                                                                           0090
C     Fortran 77 coded by F. James, 1993                                    0091
C     Converted to perform the recursion with integer arithmetic,           0092
C        K. G. Hamilton, 1996                                               0093
c     (Most of these changes are in lower case letters.)                    0094
C                                                                           0095
C     References:                                                           0096
C     M. Luscher, Computer Physics Communications  79 (1994) 100            0097
C     F. James, Computer Physics Communications 79 (1994) 111               0098
C     F. James, Computer Physics Communications 97 (1996) 357               0099
C                                                                           0100
C   LUXURY LEVELS.                                                          0101
C   ------ ------      The available luxury levels are:                     0102
C                                                                           0103
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia          0104
C           and Zaman, very long period, but fails many tests.              0105
C  level 1  (p=48): considerable improvement in quality over level 0,       0106
C           now passes the gap test, but still fails spectral test.         0107
C  level 2  (p=97): passes all known tests, but theoretically still         0108
C           defective.                                                      0109
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible             0110
C           correlations have very small chance of being observed.          0111
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.          0112
C                                                                           0113
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     0114
C!!!  Calling sequences for RANLUX:                                  ++     0115
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++     0116
C!!!                   32-bit random floating point numbers between  ++     0117
C!!!                   zero (not included) and one (also not incl.). ++     0118
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++     0119
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++     0120
C!!!               which is integer between zero and MAXLEV, or if   ++     0121
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++     0122
C!!!               should be set to zero unless restarting at a break++     0123
C!!!               point given by output of RLUXAT (see RLUXAT).     ++     0124
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++     0125
C!!!               which can be used to restart the RANLUX generator ++     0126
C!!!               at the current point by calling RLUXGO.  K1 and K2++     0127
C!!!               specify how many numbers were generated since the ++     0128
C!!!               initialization with LUX and INT.  The restarting  ++     0129
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++     0130
C!!!   A more efficient but less convenient way of restarting is by: ++     0131
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++     0132
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++     0133
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++     0134
C!!!                 32-bit integer seeds, to be used for restarting ++     0135
C!!!      ISVEC must be dimensioned 25 in the calling program        ++     0136
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     0137
      DIMENSION RVEC(LENV)                                                  0138
      DIMENSION ISEEDS(24), ISDEXT(25)                                      0139
      PARAMETER (MAXLEV=4, LXDFLT=3)                                        0140
      DIMENSION NDSKIP(0:MAXLEV)                                            0141
      DIMENSION NEXT(24)                                                    0142
      integer   iuni, icarry, maskhi, masklo                                0143
      PARAMETER (IGIGA=1000000000, JSDFLT=314159265)                        0144
      PARAMETER (ITWO24=16777216, ICONS=2147483563)                         0145
      parameter (masklo=itwo24-1)                                           0146
      SAVE NOTYET, I24, J24, icarry, ISEEDS, TWOM24, TWOM12, LUXLEV         0147
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED                 0148
      INTEGER LUXLEV                                                        0149
      LOGICAL NOTYET                                                        0150
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/      0151
      DATA I24,J24,icarry/24,10,0/                                          0152
C                               default                                     0153
C  Luxury Level   0     1     2   *3*    4                                  0154
      DATA NDSKIP/0,   24,   73,  199,  365 /                               0155
Corresponds to p=24    48    97   223   389                                 0156
C     time factor 1     2     3     6    10   on slow workstation           0157
C                 1    1.5    2     3     5   on fast mainframe             0158
c                                                                           0159
c     MASKLO is a word in which the rightmost 24 bits are set and the       0160
c     higher ones are clear.  MASKHI is MASKLO's complement: one with       0161
c     24 clear bits on the right, and all the higher ones set.              0162
c     If using 32-bit words, MASKLO=z'00FFFFFF', and MASKHI=z'FF000000'.    0163
c     Thus, ANDing with MASKLO gives a 24-bit reduced result,               0164
c     while ANDing with MASKHI indicates if there was a "carry."            0165
c                                                                           0166
      maskhi = not(masklo)                                                  0167
C                                                                           0168
C     NOTYET is .TRUE. if no initialization has been performed yet.         0169
C     Default Initialization by Multiplicative Congruential                 0170
c                                                                           0171
      IF (NOTYET) THEN                                                      0172
         NOTYET = .FALSE.                                                   0173
         JSEED = JSDFLT                                                     0174
         INSEED = JSEED                                                     0175
         WRITE(0,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED        0176
         LUXLEV = LXDFLT                                                    0177
         NSKIP = NDSKIP(LUXLEV)                                             0178
         LP = NSKIP + 24                                                    0179
         IN24 = 0                                                           0180
         KOUNT = 0                                                          0181
         MKOUNT = 0                                                         0182
         WRITE(0,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',        0183
     +       LUXLEV,'      p =',LP                                          0184
         TWOM24 = 1.                                                        0185
         DO 25 I= 1, 24                                                     0186
           TWOM24 = TWOM24 * 0.5                                            0187
           K = JSEED/53668                                                  0188
           JSEED = 40014*(JSEED-K*53668)-K*12211                            0189
           IF (JSEED .LT. 0)  JSEED = JSEED+ICONS                           0190
           ISEEDS(I) = MOD(JSEED,ITWO24)                                    0191
           NEXT(I) = I-1                                                    0192
   25    CONTINUE                                                           0193
         TWOM12 = TWOM24 * 4096.                                            0194
         NEXT(1) = 24                                                       0195
         I24 = 24                                                           0196
         J24 = 10                                                           0197
         icarry = 0                                                         0198
         if (iand(iseeds(24),maskhi).ne.0) icarry = 1                       0199
      ENDIF                                                                 0200
C                                                                           0201
C     The Generator proper: "Subtract-with-borrow",                         0202
C     as proposed by Marsaglia and Zaman,                                   0203
C     Florida State University, March, 1989                                 0204
C                                                                           0205
      DO 100 IVEC= 1, LENV                                                  0206
      iuni = iseeds(j24) - iseeds(i24) - icarry                             0207
      if (iand(iuni,maskhi).ne.0) then                                      0208
         iuni = iand(iuni,masklo)                                           0209
         icarry = 1                                                         0210
      else                                                                  0211
         icarry = 0                                                         0212
      endif                                                                 0213
      iseeds(i24) = iuni                                                    0214
      I24 = NEXT(I24)                                                       0215
      J24 = NEXT(J24)                                                       0216
      uni = float(iuni)*twom24                                              0217
c                                                                           0218
C     Small numbers (with less than 12 "significant" bits) are "padded,"    0219
C     and zero is forbidden in case someone takes a logarithm               0220
c                                                                           0221
      if (uni .lt. twom12) then                                             0222
         uni = uni + (float(iseeds(j24))*twom24)*twom24                     0223
         if (uni .le. 0.0) uni = twom24*twom24                              0224
      endif                                                                 0225
      rvec(ivec) = uni                                                      0226
c                                                                           0227
C     Skipping to luxury.  As proposed by Martin Luscher.                   0228
c                                                                           0229
      IN24 = IN24 + 1                                                       0230
      IF (IN24 .EQ. 24)  THEN                                               0231
         IN24 = 0                                                           0232
         KOUNT = KOUNT + NSKIP                                              0233
         DO 90 ISK= 1, NSKIP                                                0234
         iuni = iseeds(j24) - iseeds(i24) - icarry                          0235
         if (iand(iuni,maskhi).ne.0) then                                   0236
           iuni = iand(iuni,masklo)                                         0237
           icarry = 1                                                       0238
         else                                                               0239
           icarry = 0                                                       0240
         endif                                                              0241
         iseeds(i24) = iuni                                                 0242
         I24 = NEXT(I24)                                                    0243
         J24 = NEXT(J24)                                                    0244
   90    CONTINUE                                                           0245
      ENDIF                                                                 0246
  100 CONTINUE                                                              0247
      KOUNT = KOUNT + LENV                                                  0248
      IF (KOUNT .GE. IGIGA)  THEN                                           0249
         MKOUNT = MKOUNT + 1                                                0250
         KOUNT = KOUNT - IGIGA                                              0251
      ENDIF                                                                 0252
      RETURN                                                                0253
C                                                                           0254
C--------------------------------------------------------------------       0255
C     Entry to input integer seeds from previous run                        0256
C--------------------------------------------------------------------       0257
C                                                                           0258
      ENTRY RLUXIN(ISDEXT)                                                  0259
      NOTYET = .FALSE.                                                      0260
      TWOM24 = 1.                                                           0261
      DO 195 I= 1, 24                                                       0262
      NEXT(I) = I-1                                                         0263
  195 TWOM24 = TWOM24 * 0.5                                                 0264
      NEXT(1) = 24                                                          0265
      TWOM12 = TWOM24 * 4096.                                               0266
      WRITE(0,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'     0267
      WRITE(0,'(5X,5I12)') ISDEXT                                           0268
      DO 200 I= 1, 24                                                       0269
      iseeds(i) = isdext(i)                                                 0270
  200 CONTINUE                                                              0271
      icarry = 0                                                            0272
      if (isdext(25) .lt. 0) icarry = 1                                     0273
      ISD = IABS(ISDEXT(25))                                                0274
      I24 = MOD(ISD,100)                                                    0275
      ISD = ISD/100                                                         0276
      J24 = MOD(ISD,100)                                                    0277
      ISD = ISD/100                                                         0278
      IN24 = MOD(ISD,100)                                                   0279
      ISD = ISD/100                                                         0280
      LUXLEV = ISD                                                          0281
      IF (LUXLEV .LE. MAXLEV) THEN                                          0282
         NSKIP = NDSKIP(LUXLEV)                                             0283
         WRITE (0,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',      0284
     +                       LUXLEV                                         0285
      ELSE  IF (LUXLEV .GE. 24) THEN                                        0286
         NSKIP = LUXLEV - 24                                                0287
         WRITE (0,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV      0288
      ELSE                                                                  0289
         NSKIP = NDSKIP(MAXLEV)                                             0290
         WRITE (0,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV        0291
         LUXLEV = MAXLEV                                                    0292
      ENDIF                                                                 0293
      INSEED = -1                                                           0294
      RETURN                                                                0295
C                                                                           0296
C--------------------------------------------------------------------       0297
C     Entry to ouput seeds as integers                                      0298
C--------------------------------------------------------------------       0299
C                                                                           0300
      ENTRY RLUXUT(ISDEXT)                                                  0301
      DO 300 I= 1, 24                                                       0302
        isdext(i) = iseeds(i)                                               0303
  300 CONTINUE                                                              0304
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV              0305
      if (icarry .ne. 0.)  isdext(25) = -isdext(25)                         0306
      RETURN                                                                0307
C                                                                           0308
C--------------------------------------------------------------------       0309
C     Entry to output the "convenient" restart point                        0310
C--------------------------------------------------------------------       0311
C                                                                           0312
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)                                        0313
      LOUT = LUXLEV                                                         0314
      INOUT = INSEED                                                        0315
      K1 = KOUNT                                                            0316
      K2 = MKOUNT                                                           0317
      RETURN                                                                0318
C                                                                           0319
C--------------------------------------------------------------------       0320
C     Entry to initialize from one or three integers                        0321
C--------------------------------------------------------------------       0322
C                                                                           0323
      ENTRY RLUXGO(LUX,INS,K1,K2)                                           0324
c                                                                           0325
c     LUX either selects a luxury level, or specifies a skipping            0326
c     number (p) directly.                                                  0327
c                                                                           0328
      IF (LUX .LT. 0) THEN                                                  0329
         LUXLEV = LXDFLT                                                    0330
      ELSE IF (LUX .LE. MAXLEV) THEN                                        0331
         LUXLEV = LUX                                                       0332
      ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN                         0333
         LUXLEV = MAXLEV                                                    0334
         WRITE (0,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX           0335
      ELSE                                                                  0336
         LUXLEV = LUX                                                       0337
         DO 310 ILX= 0, MAXLEV                                              0338
           IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX                       0339
  310    CONTINUE                                                           0340
      ENDIF                                                                 0341
      IF (LUXLEV .LE. MAXLEV)  THEN                                         0342
         NSKIP = NDSKIP(LUXLEV)                                             0343
         WRITE(0,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',     0344
     +        LUXLEV,'     P=', NSKIP+24                                    0345
      ELSE                                                                  0346
          NSKIP = LUXLEV - 24                                               0347
          WRITE (0,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV     0348
      ENDIF                                                                 0349
      IN24 = 0                                                              0350
c                                                                           0351
c     INS can supply a seed value.  If INS=0, default seed is used.         0352
c                                                                           0353
      IF (INS .LT. 0)  WRITE (0,'(A)')                                      0354
     +   ' Illegal initialization by RLUXGO, negative input seed'           0355
      IF (INS .GT. 0)  THEN                                                 0356
        JSEED = INS                                                         0357
        WRITE(0,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',     0358
     +      JSEED, K1,K2                                                    0359
      ELSE                                                                  0360
        JSEED = JSDFLT                                                      0361
        WRITE(0,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'     0362
      ENDIF                                                                 0363
      INSEED = JSEED                                                        0364
      NOTYET = .FALSE.                                                      0365
      TWOM24 = 1.                                                           0366
      DO 325 I= 1, 24                                                       0367
        TWOM24 = TWOM24 * 0.5                                               0368
        K = JSEED/53668                                                     0369
        JSEED = 40014*(JSEED-K*53668)-K*12211                               0370
        IF (JSEED .LT. 0)  JSEED = JSEED+ICONS                              0371
        ISEEDS(I) = MOD(JSEED,ITWO24)                                       0372
        NEXT(I) = I-1                                                       0373
  325 CONTINUE                                                              0374
      TWOM12 = TWOM24 * 4096.                                               0375
      NEXT(1) = 24                                                          0376
      I24 = 24                                                              0377
      J24 = 10                                                              0378
      icarry = 0                                                            0379
      if (iseeds(24) .eq. 0) icarry = 1                                     0380
c                                                                           0381
C     If restarting at a break point, skip K1 + IGIGA*K2                    0382
C     Note that this is the number of numbers delivered to                  0383
C     the user PLUS the number skipped (if luxury .GT. 0).                  0384
c                                                                           0385
      KOUNT = K1                                                            0386
      MKOUNT = K2                                                           0387
      IF (K1+K2 .NE. 0)  THEN                                               0388
        maskhi = not(masklo)                                                0389
        DO 500 IOUTER=1,K2+1                                                0390
          INNER = IGIGA                                                     0391
          IF (IOUTER .EQ. K2+1)  INNER = K1                                 0392
          DO 450 ISK= 1, INNER                                              0393
            iuni = iseeds(j24) - iseeds(i24) - icarry                       0394
            if (iand(iuni,maskhi).ne.0) then                                0395
              iuni = iand(iuni,masklo)                                      0396
              icarry = 1                                                    0397
            else                                                            0398
              icarry = 0                                                    0399
            endif                                                           0400
            iseeds(i24) = iuni                                              0401
            I24 = NEXT(I24)                                                 0402
            J24 = NEXT(J24)                                                 0403
  450     CONTINUE                                                          0404
  500   CONTINUE                                                            0405
C         Get the right value of IN24 by direct calculation                 0406
        IN24 = MOD(KOUNT, NSKIP+24)                                         0407
        IF (MKOUNT .GT. 0)  THEN                                            0408
           IZIP = MOD(IGIGA, NSKIP+24)                                      0409
           IZIP2 = MKOUNT*IZIP + IN24                                       0410
           IN24 = MOD(IZIP2, NSKIP+24)                                      0411
        ENDIF                                                               0412
C       Now IN24 had better be between zero and 23 inclusive                0413
        IF (IN24 .GT. 23) THEN                                              0414
           WRITE (0,'(A/A,3I11,A,I5)')                                      0415
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,         0416
     +     K1, K2, ' cannot occur at luxury level', LUXLEV                  0417
           IN24 = 0                                                         0418
        ENDIF                                                               0419
      ENDIF                                                                 0420
      RETURN                                                                0421
      END                                                                   0422

