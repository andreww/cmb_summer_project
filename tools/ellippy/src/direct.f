	subroutine direct()
	real t1(6),t2(6),t3(6)
        character*8 phcod
        open(15,file='ELCOR.dat',action='READ')
        open(8,file='elcordir.tbl',access='direct',
     ^       recl=80,form='formatted')
        nr=1
10      read(15,*,end=11) phcod,n,d1,d2
        write(6,*) phcod,nr,n
        write(8,rec=nr,fmt='(a8,i10,2f10.0)') phcod,n,d1,d2
        nr = nr+1
        do j=1,n
         read(15,*) del
         read(15,*) (t1(m),m=1,6)
         read(15,*) (t2(m),m=1,6)
         read(15,*) (t3(m),m=1,6)
         write(8,rec=nr,fmt='(f10.0)') del
         nr = nr+1
         write(8,rec=nr,fmt='(6f10.4)') (t1(m),m=1,6)
         nr = nr+1
         write(8,rec=nr,fmt='(6f10.4)') (t2(m),m=1,6)
         nr = nr+1
         write(8,rec=nr,fmt='(6f10.4)') (t3(m),m=1,6)
         nr = nr+1
        enddo
        go to 10
11	continue
        close(15)
        close(8)
        end
