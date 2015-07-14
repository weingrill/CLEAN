	character DFILE*80,EFILE*80
	write (*,*)'what is the name of the input data file (time,data)'
	read (*,1) DFILE
 1	 format (A80)

	write (*,*)'what is the period to be removed from the data?'
	read (*,*) period
	write (*,*)'what is the name of the output data file (phase,data)'
	read (*,1) EFILE
C  open the file & read in the header and data unformatted (time,data)
        OPEN(UNIT=10,FILE=DFILE,STATUS='OLD',
     : ACCESS='SEQUENTIAL', FORM='FORMATTED')
CCARRIAGECONTROL='LIST')
        OPEN(UNIT=11,FILE=EFILE,STATUS='unknown',
     : ACCESS='SEQUENTIAL', FORM='FORMATTED')
CCARRIAGECONTROL='LIST')
	write (*,*)'skipping one line...'
	read (10,*)
 99	read (10,*,end=999) ut,data
	phase = ut/period - aint ( ut/period )
	write (11,2) phase,data
	write (11,2) phase + 1.000 , data
 2	 format (1x,F8.6,3x,F8.5)
 	go to 99
 999    continue
        write (*,*) 'finished...'
	end
