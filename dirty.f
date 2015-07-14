       PROGRAM DIRTY
       IMPLICIT NONE
       REAL      t,x,x1,x2
       CHARACTER dirtyfile*80,header*80,amplitudefile*80
       write(*,*)'enter file name of dirty spectrum?'
       read (*,'(A)') dirtyfile
       OPEN (10,file=dirtyfile,status='old')
       write(*,*)'enter file name to put amplitudes in?'
       read (*,'(A)') amplitudefile
       OPEN (11,file=amplitudefile,status='unknown')
	read (10,1)header
       write (11,1) header
 1      FORMAT (A80)
 99	read (10,*,end=999) t,x1,x2
	x=sqrt(x1**2+x2**2)
	write (11,*) t,x
 	go to 99
 999    continue
        write (*,*) 'finished...'
       stop
	end
