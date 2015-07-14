        PROGRAM DFOURT
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  NOTE: This code was written for clarity, not for efficiency.
C        Some optimization may be desired for certain applications.
C----------------------------------------------------------------------------
C  Produces the dirty spectrum and the spectral window for a time series
C  contained in XFILE.  The time average and the data mean are automatically
C  removed (see paper).  The frequency sampling is determined by the user.
C              / freq. resolution = 1/T   ; T is the overall time interval
C   defaults: <  Max. frequency   = 1/2dt ; dt is the smallest time spacing
C              \ Points-per-beam  = 4     ; Freq. spacing, dF = 0.25*(1/T)
C  The dirty spectrum is written to DFILE, and the spectral window to WFILE.
C  NOTE: the filenames are limited to 14 characters.  This is so they can
C        be propagated to subsequent programs through the spec. file headers
C----------------------------------------------------------------------------
C  CALLS:     RDDATA       to read the time series data  \_ in READWRITE.FOR
C             WRSPEC       to write spectra              /
C             DFOUR        performs the discrete FT
C----------------------------------------------------------------------------
C        PARAMETER MAXN= 5000           ! maximum # data samples
C        PARAMETER MAXM= 50000         ! maximum index for FT array
        INTEGER   MAXN, MAXM
        PARAMETER ( MAXN = 5000 )
        PARAMETER ( MAXM = 50000 )

 
C  Declare arrays and some variables
        COMPLEX      D(0:MAXM),W(0:2*MAXM),DFOUR
        REAL         T(MAXN),X(MAXN),ONES(MAXN),F(0:2*MAXM)
        CHARACTER    XFILE*25,DFILE*25,WFILE*25,HEADER*80
        DATA ONES/MAXN*1.0/           ! fill the ONES array with ones
 
C  Get execution parameters and filenames
        WRITE(*,*) 'INPUT FILES & PARAMETERS:'
        WRITE(*,*) '(''$Time series data file: '')'
        READ (5,'(A)') XFILE
        WRITE(*,*) '(''$[0=default]  Frequency resolution: '')'
        READ (5,*) FRES
        WRITE(*,*) '(''$[0=default]  Maximum frequency (for D): '')'
        READ (5,*) FMAX
        WRITE(*,*) '(''$[0=default]  Points-per-beam: '')'
        READ (5,*) PPB
        WRITE(*,*) 'OUTPUT FILES:'
        WRITE(*,*) '(''$Dirty spectrum file: '')'
        READ (5,'(A)') DFILE
        WRITE(*,*)  '(''$Spectral window file: '')'
        READ (5,'(A)') WFILE
 
C  Read the times T(1:N) and the data X(1:N) from XFILE
        N= MAXN                               ! defines max. N for RDDATA
        CALL RDDATA(N,T,X,XFILE,HEADER)
        WRITE(*,*) HEADER
	WRITE(*,*) 'total points is ',n
 
C  Find the mean T,X
        TMEAN= 0.
        XMEAN= 0.
        DO I=1,N
           TMEAN= TMEAN+ T(I)
           XMEAN= XMEAN+ X(I)
        ENDDO
        TMEAN= TMEAN/N
        XMEAN= XMEAN/N
	write (*,*) tmean,xmean
 
C  subtract the time,data means from T,X
        DO I=1,N
           T(I)= T(I)-TMEAN
           X(I)= X(I)-XMEAN
		write (12,*) t(i),x(i)
        ENDDO
 
C  If parameters are supplied as zero, set default values
        SMIN= 1.E20                          !\--(larger than expected SEPs)
        DO I=2,N                             ! \
           SEP= T(I) - T(I-1)                !  > find min. time separation
           IF(SEP.LT.SMIN) SMIN= SEP         ! /
        ENDDO                                !/
        SMAX= T(N)-T(1)                      !> max. time separation
        IF(FRES.EQ.0.) FRES= 1./SMAX         ! frequency resolution
        IF(FMAX.EQ.0.) FMAX= 1./(2.*SMIN)    ! max. frequency
        IF(PPB.EQ.0.) PPB= 4.                ! points-per-beam
C    ---confirm parameter values
        WRITE(*,*) 'FRES,FMAX,PPB=',FRES,FMAX,PPB
 
C  set up the frequency array F(0:M)
        dF= FRES/PPB                   ! frequency increment
        M= INT(FMAX/dF)                ! maximum freq. element (for D)
        DO J=0,2*M
           F(J)= dF*J
        ENDDO
 
C  calculate the dirty spectrum D(0:M) and the spec. window W(0:2M)
        WRITE(*,*) 'Computing the dirty spectrum...'
        DO J= 0,M
           D(J)= DFOUR(F(J),N,T,X)
        ENDDO
        WRITE(*,*) 'Computing the spectral window...'
        DO J= 0,2*M
           W(J)= DFOUR(F(J),N,T,ONES)
        ENDDO
 
C  create the headers & write to the output files
C    ---dirty spectrum (Header includes N, Tsub, Xsub for later use)'
        WRITE(HEADER,'(A14,'',df='',A14,'',N='',I5,'',Tsub='',
     &   1P,E14.7,'',Dsub='',1P,E14.7)') DFILE,XFILE,N,TMEAN,XMEAN
        CALL WRSPEC(M,F,D,DFILE,HEADER)
        WRITE(*,*) HEADER
C    ---spectral window (Header includes N, Tsub, Tmean for later use)'
        WRITE(HEADER,'(A14,'',df='',A14,'',N='',I5,'',Tsub='',
     &   1P,E14.7,'',Tavg='',1P,E14.7)') WFILE,XFILE,N,TMEAN,TMEAN
        CALL WRSPEC(2*M,F,W,WFILE,HEADER)
        WRITE(*,*) HEADER
 
C  exit the program
        CALL EXIT
        END
 
 
 
        COMPLEX FUNCTION DFOUR(FREQ,N,TIME,DATA)
C----------------------------------------------------------------------------
C  Author:   J. Lehar                                         Date: 13-JUL-87
C----------------------------------------------------------------------------
C  Returns the Fourier transform of the time series specified by TIME(1:N)
C  and DATA(1:N), evaluated at the frequency FREQ.
C  The form of the Fourier transform is taken from Bracewell (1965)
C
C                     1  .--. N           -i*2*PI*F*TIME(i)
C           DFT(F) = ---  >      DATA(i) e
C                     N  `--'i=1
C
C  The DFT is normalized to have the data mean at FREQ=0.
C----------------------------------------------------------------------------
        COMPLEX    SUM
        REAL       DATA(N),TIME(N)
 
C  initialize some variables
        SUM= (0.,0.)                    ! accumulation variable
        TWOPI= 8.*ATAN2(1.,1.)          ! calculate 2*pi
 
C  Evaluate FT at F...
        DO I=1,N
           PHASE= -TWOPI*FREQ*TIME(I)
           SUM= SUM + DATA(I)*CMPLX(COS(PHASE),SIN(PHASE))
        ENDDO
 
C  return with FT properly normalized
        DFOUR= SUM/N           ! ensures correct normalization
        RETURN
        END
 
 
