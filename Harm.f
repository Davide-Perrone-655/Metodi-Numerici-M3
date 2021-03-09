	program Harmonic
	
	!ricorda di avere L*eta costante nelle varie scale, corrispondono ad w e beta
	integer, parameter :: L=10000
	integer, parameter :: meas=1000
	integer, parameter :: term=100
	real :: lattice(L)	
	
	integer :: seed
	real :: eta, delta, bhw
	
	integer :: bcp(L), bcm(L)
	integer :: warm, i, j, s
	
	real :: misura1, temp
	
	common/parameters/eta, delta
	common/ranseed/seed
	
	!definisco tutti i parametri iniziali
	seed=124
	bhw=3.
	eta=bhw/L
	!delta=0.5
	delta=sqrt(eta)
	warm=meas/10						!quanti giri a vuoto faccio
	
	!stampo i parametri
	write(*,*) 'parametri scelti'
	
	write(*,*) 'seed: ', seed
	write(*,*) 'lunghezza: ', L
	write(*,*) 'beta hbar omega', bhw
	write(*,*) 'eta: ', eta
	write(*,*) 'delta: ', delta
	
	write(*,*) 'misure: ', meas
	write(*,*) 'warm-up: ', warm
	write(*,*) 'Random step: NO '
	
	
	
	
	call boundary_geom(bcp, bcm, L)
	call initialize(lattice,L)
	call warmup(lattice, bcp, bcm, L, warm)
	
	call saveparam(meas)
	
	
	!open(10, file='storia_prova.dat')
	!open(11, file='misure_prova.dat')
	
	do i=1, meas
		do j=1, term
			call chainstep(lattice, bcp, bcm, L)
		end do
		!write(10,*) lattice
		!call misurazione(lattice, L, misura1)
		!write(11,*) misura1
	end do
	
	
	
	!call savepath(lattice, L)
	
	
	
	CONTAINS
	
	
c=========================================================
	!INIZIALIZZAZIONI
C=========================================================	

	!Condizioni al bordo
	subroutine boundary_geom(bcp, bcm, len_row)
	integer :: len_row, bcp(len_row), bcm(len_row)
	integer :: i
	
	do i=1,len_row
		bcm(i)=i-1
		bcp(i)=i+1
	end do
	bcm(1)=len_row
	bcp(len_row)=1
	
	!write(*,*) bcm
	!write(*,*) bcp
	end subroutine boundary_geom
	
	
	
	

	!inizializza la matrice di spin in modo random
	subroutine initialize(vector, len_row)
	common /ranseed /seed
	integer :: len_row, seed
	real :: vector(len_row)
	
	integer :: i
	real :: num
	
	do i=1,len_row
		num= 2*ran2(seed) - 1			!inizializzo con un numero random tra -1 e 1
		vector(i)=num
	end do
	!write(*,*) vector
	end subroutine initialize
	
	
	
	
	!funzione per riscaldare il path integral (faccio giri a vuoto)
	subroutine warmup(vector, bcp, bcm, len_row, w_len)			!gli passo le cose del chainstep + w_len: numero di step da fare 
	common/ranseed/seed
	common/parameters/eta, delta
	integer :: bcp(len_row), bcm(len_row)
	integer :: seed, len_row,  w_len
	real :: eta, delta
	real :: vector(len_row)
	integer :: i
	
	do i=1, w_len
		call chainstep(vector, bcp, bcm, len_row)
	end do
	end subroutine warmup
	
	
	subroutine saveparam(len_row)
	common/parameters/eta, delta
	integer :: len_row
	real :: delta, eta
	open(115, file='shape_10.dat')
	write(115,*) len_row, delta, eta
	close(115)
	end subroutine saveparam
	
	
	
C==========================================================	
	!PASSO BASE METROPOLIS
C==========================================================

	!passo non random
	subroutine chainstep(vec, bcp, bcm, len_row)
	common/parameters/eta, delta
	common/ranseed/seed
	
	integer :: seed
	real :: eta, delta
	
	real :: vec(len_row)
	integer ::  bcp(len_row), bcm(len_row)
	integer :: len_row
	
	real :: z, f1, yp, y, r
	real :: c1, c2
	integer :: i0, im, ip, i
	
	
	c2= ((eta/2.) + 1./eta)
	c1=1./eta
	
	
	do i=1,len_row
		y=vec(i)								
		!chiamo le funzioni per assicurare la corretta geometria sui bordi
		im=bcm(i)
		ip=bcp(i)


		yp= y + 2.*delta*(0.5 - ran2(seed))				!cerco un nuovo valore, uniformemente in un intervallo di raggio delta
		f1=vec(im) + vec(ip)
		
		r=  c1*yp*f1 - c2*yp**2
		r= r - c1*y*f1 + c2*y**2 
		
		
		!calcolo la probabilità di accettare il valore
			
		z=log(ran2(seed))							!elemento per fare il confronto

		if(z.lt.r) vec(i)=yp			!test metropolis
		
	end do
	
	end subroutine chainstep
	
	
	
	!passo random
	subroutine ranstep(vec, bcp, bcm, len_row)
	common/parameters/eta, delta
	common/ranseed/seed

	
	integer :: seed
	real :: eta, delta
	
	real :: vec(len_row)
	integer ::  bcp(len_row), bcm(len_row)
	integer :: len_row
	
	real :: z, f0, f1, yp, y, r
	real :: c1, c2
	integer :: i0, im, ip, i
	
	
	c1= ((eta/2.) + 1./eta)
	c2=1./eta
	
	do i=1,len_row
		
		i0=ceiling(ran2(seed)*(len_row))					!prendo un elemento random
		y=vec(i0)								
		
		z=log(ran2(seed))							!elemento per fare il confronto
		
		
		!chiamo le funzioni per assicurare la corretta geometria sui bordi
		im=bcm(i0)
		ip=bcp(i0)
		
		yp= y + delta*(1.0 - 2.0*ran2(seed))
			f1=(vec(im)+vec(ip))
		f0= - c1*(yp**2) + c2*yp*f1
		r=f0 + c1*(y**2) - c2*y*f1
		
		
			!calcolo la probabilità di accettare il valore
		if(z.lt.r) vec(i0)=yp			!accetto se z è più piccolo di r 
	end do
	
	end subroutine ranstep
	
C==============================================================================
C	Funzioni varie
C==============================================================================

	!funzione per salvare il cammino (non testata)
	subroutine savepath(vector, len_row)
	common/parameters/eta, delta
	real :: vector(len_row), eta, delta
	integer :: len_row
	integer :: i
	
	open(133, file='Path_10000')
	!write(133,*) '#len', 'eta', 'delta'
	!write(133,*) '#', len_row, eta, delta
	
	do i=1, len_row
		write(133,*) vector(i)
	end do
	close(133)
	end subroutine savepath
	




	subroutine misurazione(vector, len_row, res)
	real :: vector(len_row), res
	integer :: len_row
	
	real :: obs
	integer :: i
	
	obs=0.0
	do i=1, len_row
		obs=obs + vector(i)**2
	end do
	
	res=obs/float(len_row)
	
	end subroutine misurazione
	
	
	
	
	
	
c==============================================================
C	Random number generator: ran2 di numerical recipes
c==============================================================
	FUNCTION ran2(idum)
	INTEGER :: idum
	REAL :: ran2
	INTEGER, PARAMETER :: IM1=2147483563, IM2=2147483399, IR1=12211
	INTEGER, PARAMETER :: IR2=3791, NTAB=32,  IMM1=IM1-1
	INTEGER,PARAMETER :: IA1=40014,IA2=40692,IQ1=53668,IQ2=52774
	real, parameter :: AM=1./IM1, NDIV=1+IMM1/NTAB, EPS=1.2e-7, RNMX=1.-EPS
	
	INTEGER :: idum2,j,k,iv(NTAB),iy 
	SAVE iv,iy,idum2
	DATA idum2/123456789/, iv/NTAB*0/, iy/0/
	if (idum.le.0) then
		idum=max(-idum,1)
		idum2=idum
		do j=NTAB+8,1,-1
			k=idum/IQ1
			idum=IA1*(idum-k*IQ1)-k*IR1
			if (idum.lt.0) idum=idum+IM1
			if (j.le.NTAB) iv(j)=idum
		end do
		iy=iv(1)
	endif
	k=idum/IQ1
	idum=IA1*(idum-k*IQ1)-k*IR1
	if (idum.lt.0) idum=idum+IM1
	k=idum2/IQ2
	idum2=IA2*(idum2-k*IQ2)-k*IR2
	if (idum2.lt.0) idum2=idum2+IM2
	j=1+iy/NDIV
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+IMM1
	ran2=min(AM*iy,RNMX)
	return
	end function
	
	
	end program harmonic
