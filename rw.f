        SUBROUTINE RDDATA(N,TIME,DATA,DFILE,HEADER)
C----------------------------------------------------------------------------
C  Author:  J. Lehar                                         Date:  13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  Source file:  READWRITE.FOR
C----------------------------------------------------------------------------
C  Reads in the times, data values for a time series from DFILE
C  In addition, a descriptive header is read in.  N (input) specifies the
C  maximum number of samples allowed and returns (output) the number of
C  samples found.  OPEN uses UNIT=1
C     ARGUMENTS:
C         N          (input) #samples allowed, (output) #samples found
C         TIME       (output) array for time samples
C         DATA       (output) array for data samples
C         DFILE      (input) data file name
C         HEADER     (output) descriptive file header
C----------------------------------------------------------------------------
        REAL TIME(N),DATA(N)
        CHARACTER DFILE*(*),HEADER*80
 
C  open the file & read in the header
        OPEN(UNIT=1,FILE=DFILE,STATUS='OLD',CARRIAGECONTROL='LIST')
        READ(1,'(A80)') HEADER
 
C  read in the time series, each record has (Time,Data)
        I=1
C       step through the file until END-OF-FILE
10         READ(1,*,END=11) TIME(I),DATA(I)
           I= I+1
           IF(I.GT.N)THEN
              TYPE *,'RDDATA: data exceeds array size, truncating at N=',N
              GOTO 11
           ENDIF
           GOTO 10
11      CONTINUE
 
C  close the file & return the number of samples found
        CLOSE(UNIT=1)
        N= I-1
        RETURN
        END
 
 
        SUBROUTINE WRDATA(N,TIME,DATA,DFILE,HEADER)
C----------------------------------------------------------------------------
C  Author:  J. Lehar                                         Date:  13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  Source file:  READWRITE.FOR
C----------------------------------------------------------------------------
C  Writes the times, data values for a time series from DFILE
C  In addition, a descriptive header is written. OPEN uses UNIT=1
C     ARGUMENTS:
C         N          (input) #samples
C         TIME       (input) array for time samples
C         DATA       (input) array for data samples
C         DFILE      (input) data file name
C         HEADER     (input) descriptive file header
C----------------------------------------------------------------------------
        REAL TIME(N),DATA(N)
        CHARACTER DFILE*(*),HEADER*80
 
C  open the file & write the header
        OPEN(UNIT=1,FILE=DFILE,STATUS='NEW',CARRIAGECONTROL='LIST')
        WRITE(1,'(A80)') HEADER
 
C  write the time series
        DO I=1,N
           WRITE(1,'(2(1X,1P,E14.7))') TIME(I),DATA(I)
        ENDDO
 
C  close the file & return
        CLOSE(UNIT=1)
        RETURN
        END
 
 
        SUBROUTINE RDSPEC(M,FREQ,SPEC,SFILE,HEADER)
C----------------------------------------------------------------------------
C  Author:  J. Lehar                                         Date:  13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  Source file:  READWRITE.FOR
C----------------------------------------------------------------------------
C  Reads in the frequency and spectral values for a spectrum to SFILE.
C  In addition, a descriptive header is read in.  M (input) specifies the
C  maximum frequency element allowed and returns (output) the maximum
C  frequency element found.  OPEN uses UNIT=1
C     ARGUMENTS:
C         M       (input) max. freq element allowed, (output) max. elt. found
C         FREQ    (output) array for frequency samples
C         SPEC    (output) array for spectrum samples
C         SFILE   (input) spectrum file name
C         HEADER  (output) descriptive file header
C----------------------------------------------------------------------------
        COMPLEX SPEC(0:M)
        REAL FREQ(0:M)
        CHARACTER SFILE*(*),HEADER*80
 
C  open the file & read in the header
        OPEN(UNIT=1,FILE=SFILE,STATUS='OLD',CARRIAGECONTROL='LIST')
        READ(1,'(A80)') HEADER
 
C  read in the spectrum, each record has (Freq,Real.comp,Imag.comp)
        I=0
C       step through the file until END-OF-FILE
10         READ(1,'(3(1X,1P,E14.7))',END=11) FREQ(I),SPEC(I)
           I= I+1
           IF(I.GT.M)THEN
              TYPE *,'RDSPEC: data exceeds array size, truncating at M=',M
              GOTO 11
           ENDIF
           GOTO 10
11      CONTINUE
 
C  close the file & return the max. frequency element found
        CLOSE(UNIT=1)
        M= I-1
        RETURN
        END
 
 
        SUBROUTINE WRSPEC(M,FREQ,SPEC,SFILE,HEADER)
C----------------------------------------------------------------------------
C  Author:  J. Lehar                                         Date:  13-JUL-87
C  Written for VAX FORTRAN 77   (May require some changes for other machines)
C  Source file:  READWRITE.FOR
C----------------------------------------------------------------------------
C  Writes the frequency and spectral values for a spectrum to SFILE
C  In addition, a descriptive header is written. OPEN uses UNIT=1
C     ARGUMENTS:
C         M       (input) max. frequency element
C         FREQ    (input) array for frequency samples
C         SPEC    (input) array for spectrum samples
C         SFILE   (input) spectrum file name
C         HEADER  (input) descriptive file header
C----------------------------------------------------------------------------
        COMPLEX SPEC(0:M)
        REAL FREQ(0:M)
        CHARACTER SFILE*(*),HEADER*80
 
C  open the file & read in the header
        OPEN(UNIT=1,FILE=SFILE,STATUS='NEW',CARRIAGECONTROL='LIST')
        WRITE(1,'(A80)') HEADER
 
C  write out the spectrum, the spec. file lists FREQ,REAL(SPEC),IMAG(SPEC)
        DO I=0,M
           WRITE(1,'(3(1X,1P,E14.7))') FREQ(I),SPEC(I)
        ENDDO
 
C  close the file & return
        CLOSE(UNIT=1)
        RETURN
        END
 
