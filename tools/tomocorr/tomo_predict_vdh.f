c===================================================================
c	EJG Jan 2004
c	use VDH 1999 tomographic model (20 layers)
c===================================================================

        character*80 infil, model1d, model3d
        character junk*1
        dimension tomo(89,180,20),top(20),bot(20)
	dimension v1d(10000)
	dtr=0.017453292
        rtd=57.29577951
        rad = 6371.0
	sumdt=0.

c . . . collect the standard input
c	"model1d" contains the reference structure
	read(*,"(a80)")model1d
c	"model3d" contains the 3D heterogeneity model
	read(*,"(a80)")model3d
c	"infil" contains the taup output of path's XYZ 
	read(*,"(a80)")infil
c       read in the depth range for calculating the tomog prediction:
        read(*,*)ztop,zbot

c . . . load in model1D: reference structure 
c	note: make sure the 3D structure is w.r.t. the 1D model!
	open(11,file=model1d)
	i=1
10	read(11,*,end=12)z,v1d(i)
	i=i+1
	goto10
12	n_v1d=i-1
	close(11)
	
c . . . load in model3D: tomography model 
	open(11,file=model3d)
	do 18 layer=1,20
		read(11,*)top(layer),bot(layer)
  		do 18 l=89,1,-1
			do 18 k=1,180
18			read(11,*)tomo(l,k,layer)
	close(11)


c . . . open file to print pathstats:
	open(23,file="pathstats")
	write(23,*)'    delta   z        ARC     v    residual  lat    lon      dt_pred   sumdt'
        write(23,*)'    (deg)  (km)     (km)   (km/s)   (%)    (deg)  (deg)      (sec)    (sec)'

c . . . open file containing path of lats/lons:
	open(15,file=infil)
	read(15,"(a1)")junk
        read(15,*)delta,R,plat,plon
        xlast = R
        ylast = 0.

20      read(15,*,end=30)delta,R,plat,plon
	if(R.gt.6370.)R=6370.
	DEL = delta*dtr
	x = R * cos(DEL)
	y = R * sin(DEL)
        ARC=sqrt((x-xlast)**2 + (y-ylast)**2)
	z = rad - R
        v = v1d(int(z))
        if(z.ge.ztop .and. z.le.zbot)then
	   call get_dt(plat,plon,z,v,ARC,dt_pred,residual,top,bot,tomo)
        else
           dt_pred = 0.0
	   residual = 0.0
        endif
 	sumdt=sumdt+dt_pred
        write(23,25)delta,z,ARC,v,residual,plat,plon,dt_pred,sumdt
25	format(7f8.2,5f10.4)
        xlast=x
        ylast=y
	goto 20
30	close(15)
	close(23)
 	write(*,*)sumdt
	end 

c===================================================================
	subroutine get_dt(lat1,lon1,z,VPREM,sddp,dt,resid,top,bot,tomo)
c===================================================================
        REAL lat1,lon1
        dimension tomo(89,180,20),top(20),bot(20)

c	get VDH model depth index 
	do 300 izi=1,20
		if(z.ge.top(izi).and.z.lt.bot(izi))then
			layer = izi
			goto 301
		endif
300	continue

c	get VDH grid pt coords surrounding lat,lon
301	lathi	= int(lat1/2.)*2
	if(lat1.gt.0.) then
		latlo=lathi
		lathi=lathi + 2
	else
		latlo = lathi - 2
	endif
	if(latlo.lt.-88) latlo=-88
	if(lathi.gt.88) lathi=88
	lonhi = int(lon1/2.)*2
	if(lon1.gt.0 ) then
		lonlo=lonhi
		lonhi=lonhi + 2
	else
		lonlo = lonhi - 2
	endif
	lonlo = lonhi - 2
	ilatlo = 45 + latlo/2
	ilathi = ilatlo + 1
	ilonlo = 90 + lonlo/2
	ilonhi = ilonlo + 1

c	get resids @ corners, then @ lat1,lon1 (all in %)
	fact = ( lat1 - float(latlo) ) / 2.
	vlatlo=tomo(ilatlo,ilonlo,layer)*(1-fact) + fact*tomo(ilathi,ilonlo,layer)
	vlathi=tomo(ilatlo,ilonhi,layer)*(1-fact) + fact*tomo(ilathi,ilonhi,layer)
	fact = ( lon1 - float(lonlo) ) / 2.
	resid = vlatlo + fact * ( vlathi - vlatlo )

c	get time delay assoc. with resid
	Tprem = sddp/VPREM
	Tanom = sddp / ( VPREM * ( 1 + resid/100. ) )
	dt = Tanom - Tprem
	return
	end
