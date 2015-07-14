        PROGRAM CLEAN
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  NOTE: This code was written for clarity, not for efficiency.
C        Some optimization may be desired for certain applications.
C----------------------------------------------------------------------------
C  Deconvolves the spectral window in WFILE from the dirty spectrum in DFILE,
C  by using an iterative, one-dimensional, complex, Hoegbom CLEAN algorithm.
C  The resulting clean spectrum is in SFILE, the residual and CLEAN component
C  spectra are found in RFILE & CFILE, respectively.
C     The user determines the gain and the number of cleans, the latter
C  either explicitly or by specifying a CLEAN level.  (when the maximum in
C  the residual spectrum is less than the CLEAN level, cleaning stops.)
C    Since all spectra from real data are Hermitian, we define only the non-
C  negative frequencies.  The negative frequency elements are recovered by
C  the use of the function CVAL, which returns the complex conjugate for
C  negative frequencies, and zero for frequencies outside the defined range.
C  NOTE: the filenames are limited to 14 characters.  This is so they can
C        be propagated to subsequent programs through the spec. file headers
C----------------------------------------------------------------------------
C  CALLS:     RDSPEC               to read spectra   \_ in READWRITE.FOR
C             WRSPEC               to write spectra  /
C             CLEAN1               performs 1 iteration of CLEAN
C             FITBEAM              to fit a beam to the spectral window
C             RESTORE              to create the clean spectrum
C----------------------------------------------------------------------------
        INTEGER   MAXM 
        PARAMETER ( MAXM = 40000 )
C        PARAMETER MAXM= 40000      ! maximum index for spectral arrays
 
C  Declare arrays and some variables
        COMPLEX   W(0:2*MAXM),B(0:MAXM),R(0:MAXM),C(0:MAXM),S(0:MAXM)
        REAL      F(0:2*MAXM)
        CHARACTER WFILE*25,DFILE*25,RFILE*25,CFILE*25,SFILE*25,HEADER*80
 
C  get execution parameters & filenames
        WRITE(*,*) 'INPUT FILES:'
        WRITE(*,*)  '(''$Dirty spectrum file: '')'
        READ (5,'(A)') DFILE
        WRITE(*,*)  '(''$Spectral window file: '')'
        READ (5,'(A)') WFILE
        WRITE(*,*)  'INPUT PARAMETERS:'
        WRITE(*,*)  '(''$Gain: '')'
        READ (5,*) GAIN
        WRITE(*,*) '(''$[ >0 = #CLEANs, <0 = CLEAN level ] Cleans: '')'
        READ (5,*) CLEANS
        WRITE(*,*)  'OUTPUT FILES:'
        WRITE(*,*)  '(''$Residual spectrum file: '')'
        READ (5,'(A)') RFILE
        WRITE(*,*)  '(''$CLEAN components file: '')'
        READ (5,'(A)') CFILE
        WRITE(*,*)  '(''$Clean spectrum file: '')'
        READ (5,'(A)') SFILE
C    ---determine depth of CLEAN
        IF(CLEANS.GE.0.)THEN
           NCLEAN= INT(CLEANS)   ! clean to the specified # of cleans
           RCLEAN= 0.            ! no level specified.
        ELSE
           NCLEAN= 32000.        ! more than anyone is likely to want
           RCLEAN= ABS(CLEANS)   ! specifies the desired level
        ENDIF
 
C  read the input spectra, get time info for beam, determine freq. range.
C    ---dirty spectrum, read into the residual R(0:MS)
        MS= MAXM                         ! max. allowed MS, for RDSPEC
        CALL RDSPEC(MS,F,R,DFILE,HEADER)
        WRITE(*,*) HEADER
C    ---spectral window W(0:MW)
        MW= 2*MAXM                       ! max. allowed MW, for RDSPEC
        CALL RDSPEC(MW,F,W,WFILE,HEADER)
        WRITE(*,*) HEADER
        READ(HEADER,'(46X,E14.7,6X,E14.7)') TZERO,TMEAN  ! from WFILE
        dF= F(1)-F(0)                    ! frequency increment
        IF(2*MS.GT.MW) MS= MW/2          ! ensure that MS < MW/2
 
C  CLEAN the residual spectrum, storing the components in C(0:MS)
        WRITE(*,*) 'CLEANing up to ',NCLEAN,'iterations,'
        WRITE(*,*) ' or down to ', RCLEAN
        DO K=1,NCLEAN                          ! on until NCLEAN
           CALL CLEAN1(MW,W,MS,R,C,GAIN,RMAX)  ! execute a single CLEAN
           IF(RMAX.LT.RCLEAN) GOTO 10          ! skip out if down to RCLEAN
        ENDDO
10      NCLEAN= K-1                            ! actual number of CLEANs
        WRITE(*,*) 'NCLEAN=',NCLEAN,'           RMAX=',RMAX
 
C  Generate the beam B(0:MB) and the clean spectrum S(0:MS)
        MB= MS                          ! max. possible MB, for FITBEAM
        CALL FITBEAM(MB,B,MW,W,dF,TZERO,TMEAN)  ! create the beam
        CALL RESTORE(MS,S,C,R,MB,B)     ! restore S (convolve CxB, add R)
 
C  Write the output spectrum
C    ---residual spectrum
        WRITE(HEADER,1000) RFILE,DFILE,WFILE,GAIN,CLEANS
        CALL WRSPEC(MS,F,R,RFILE,HEADER)
        WRITE(*,*) HEADER
C    ---CLEAN components
        WRITE(HEADER,1000) CFILE,DFILE,WFILE,GAIN,CLEANS
        CALL WRSPEC(MS,F,C,CFILE,HEADER)
        WRITE(*,*) HEADER
C    ---clean spectrum
        WRITE(HEADER,1000) SFILE,DFILE,WFILE,GAIN,CLEANS
        CALL WRSPEC(MS,F,S,SFILE,HEADER)
        WRITE(*,*) HEADER
C     - format statement for spectrum file headers...
1000   FORMAT(A14,',From:',A14,', and:',A14,'; G=',F7.4,',C=',1P,E11.4)
 
C exit the program
        CALL EXIT
        END
 
 
        SUBROUTINE CLEAN1(MW,WIN,MS,RES,CMP,GAIN,RMAX)
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C----------------------------------------------------------------------------
C  Performs one iteration of the RLD complex, 1-D, CLEAN algorithm.
C  o  The component, CCPOS, responsible for the max., at L, in the residual
C     spectrum RES, is estimated, taking into account the interference with
C     its conjugate, CCNEG at -L.   ( through RES(L)= CCPOS +CCNEG*WIN(2L) )
C  o  The spectral window WIN is matched to these components, multiplied by
C     a GAIN factor and subtracted from the residual spectrum.
C  o  GAIN*CCPOS is added to the CMP spectrum.
C  o  returns RMAX, the max. in the residual (absolute value)
C----------------------------------------------------------------------------
C  CALLS:     CVAL                 determines spectral values beyond (0:M)
C----------------------------------------------------------------------------
        COMPLEX WIN(0:MW),RES(0:MS),CMP(0:MS),WIN2L,CCPOS,CCNEG,CVAL
        DATA ERROR/0.0001/                     ! allowed error in WNORM
 
C  find the location L of max. in the residual spectrum RES(0:M)
        RMAX= CABS(RES(0))
        L= 0
        DO J=1,MS
           IF(CABS(RES(J)).GT.RMAX)THEN
              RMAX= CABS(RES(J))
              L= J
           ENDIF
        ENDDO
 
C  estimate the component at L which yields RES(L) through WIN(L)
        WIN2L= WIN(2*L)                       ! (L,-L) interference
        WNORM= 1.0-CABS(WIN2L)**2
        IF(WNORM.LT.ERROR) CCPOS= 0.5*RES(L)  ! prevent singularities
        IF(WNORM.GE.ERROR) CCPOS= (RES(L)-WIN2L*CVAL(RES,-L,MS))/WNORM
 
C  Remove effects of both +L and -L components from RES, add CCPOS to CMP
        CCPOS= GAIN*CCPOS                     ! remove GAIN*CCPOS
        CCNEG= CONJG(CCPOS)                   ! conjugate at -L
        DO J=0,MS
           RES(J)= RES(J) -CCPOS*CVAL(WIN,J-L,MW) -CCNEG*WIN(J+L)
        ENDDO      !              (from L )          (from -L )
        CMP(L)= CMP(L) + CCPOS                ! add it to CMP
 
C  return to the calling program
        RETURN
        END
 
 
        SUBROUTINE FITBEAM(MB,BEAM,MW,WIND,dF,TZERO,TMEAN)
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C----------------------------------------------------------------------------
C  Fits a clean beam, BEAM, to the primary peak of the spectral window, WIND.
C  The modulus is a Gaussian, matched to have the same half-width,half-max.
C  through linear interpolation between frequency samples.
C  The phase is a linear function of frequency, determined by the frequency
C  increment dF and the time offset from the mean time (TZERO-TMEAN).
C  NOTE:  DFOURT subtracts the time average, so that TZERO=TMEAN always.
C         therefore, beams for files from DFOURT will always be purely real.
C         If two spectra from different time series are to be added, DFOURT
C         must be modified to allow a fixed time to be subtracted from both
C         time series. Then the phases in the spectra and the FITBEAM phases
C         will be correctly constructed for combination of the spectra.
C----------------------------------------------------------------------------
        COMPLEX BEAM(0:MB),WIND(0:MW)
 
C  Find the half-width,half-max of CABS(WIND)
        HALFMAX= 0.5*CABS(WIND(0))         ! Max. of WIND at J=0
        HALFWID= 0.                        !<---(not yet found)
        DO J=1,MB                          ! step along WIND
           IF(CABS(WIND(J)).LT.HALFMAX)THEN              !\   When .LT.HALFMAX,
              WCURR= CABS(WIND(J))                       ! \  use lin. inter-
              WLAST= CABS(WIND(J-1))                     !  > polation to find
              SLOPE= 1./(WCURR-WLAST)                    ! /  the half-width
              HALFWID= FLOAT(J-1) + SLOPE*(HALFMAX-WLAST)!/
              GOTO 10                                    ! then skip out
           ENDIF
        ENDDO                              ! keep on until up to MB
10      IF(HALFWID.LE.0) WRITE(*,*)'Could not find HALF-WIDTH,HALF-MAX'
 
C  match a Gaussian of width SIGMA to the above, determine phase increment
        TWOPI= 8.*ATAN2(1.,1.)             ! calculate 2*PI to precision
        SIGMA= HALFWID/SQRT(2.*ALOG(2.))   !\ parameters for Gaussian
        CONST= 1./(2.*SIGMA*SIGMA)         !/
        dPHASE= TWOPI*dF*(TZERO-TMEAN)     ! dPHASE from time info.
 
C  Fill the BEAM array
        MB= INT(5.*SIGMA) +1               ! only fill out to 5 sigmas
        DO J= 0,MB
           X= FLOAT(J)
           GAUSS= EXP(-CONST*X*X)          ! Gaussian envelope
           PHASE= dPHASE*X                 ! phase linear in F
           BEAM(J)= GAUSS*CMPLX(COS(PHASE),SIN(PHASE))
        ENDDO
 
C  return to the calling program
        RETURN
        END
 
 
        SUBROUTINE RESTORE(MS,SP,CC,RS,MB,BM)
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C----------------------------------------------------------------------------
C  Creates the clean spectrum SP(0:MS), by convolving the CLEAN
C  components CC(0:MS) with the clean beam BM(0:MB), and finally
C  adding the residual spectrum RS(0:MS) to the result.
C     This produces a spectrum estimate which more properly represents
C  the frequency resolution (through the beam width) and which includes
C  the fluctuations in the residual spectrum.  The CLEAN spectrum will
C  also combine components within one beam-width of each other to build
C  up a signal which might appear as several smaller ones in CC.  These
C  should be combined, because frequency separations of less than one
C  beam width are not reliable.
C----------------------------------------------------------------------------
C  CALLS:     CVAL                 determines spectral values beyond (0:M)
C----------------------------------------------------------------------------
        COMPLEX SP(0:MS),CC(0:MS),RS(0:MS),BM(0:MB),CVAL
 
C  Convolve CC with BM, then add RS, to form SP
        DO J=0,MS
           SP(J)= (0.,0.)    ! reset SP(J)
           DO K=-MB,MB
              SP(J)= SP(J)+CVAL(BM,K,MB)*CVAL(CC,J-K,MS)  ! convolution
           ENDDO
           SP(J)= SP(J)+RS(J)
        ENDDO
 
C  return to the calling program
        RETURN
        END
 
 
 
 
        COMPLEX FUNCTION CVAL(SPEC,L,M)
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C----------------------------------------------------------------------------
C  Returns the "value" of the spectrum array SPEC(0:M) at the "location" L.
C     if  L  >= 0    CVAL returns SPEC(L)
C     if  L  >  0    CVAL returns the complex conjugate of SPEC(L)
C     if |L| >  M    CVAL will return (0.,0.)
C----------------------------------------------------------------------------
        COMPLEX SPEC(0:M)
 
C  set value to conjugate of absolute location if necessary
        LOC= ABS(L)               ! where necessery information is found
        IF(LOC.GT.M)THEN                    !\
           CVAL= (0.,0.)                    ! \   if LOC exceeds range
           RETURN                           ! /   return (0.,0.)
        ENDIF                               !/
        IF(L.GE.0) CVAL= SPEC(LOC)          ! return SPEC(L)
        IF(L.LT.0) CVAL= CONJG(SPEC(LOC))   ! return SPEC(-L)
        RETURN
        END
 
