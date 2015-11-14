      character*1 lst(50,60),lin(80),numb(14)
      common lst,ifr(30,5),itab(30,5),ifrs(30,5)
      data numb/'0','1','2','3','4','5','6','7','8','9',
     * '.','C','H','I'/
      itmax=5
      ist=0
      ist1=1
      open(5,file='iter')
      read(5,*)iter,nn,kk,lines
      if(iter.gt.0)go to 2
      close(5)
      open(6,file='inp')
      do 3 j=1,lines
      read(6,1)(lst(j,k),k=1,60)
 1    format(80a1)
 3    continue
      do 5 i=1,nn
      do 5 j=1,kk
      itab(i,j)=0
 5    ifr(i,j)=0
 4    read(6,*,end=6)i1,i2
      ifr(i1,i2)=1
      go to 4
 6    close(6)
      best=100000000.
      best1=best
 99   iter=iter+1
      open(5,file='iter')
      write(5,10)iter,nn,kk,lines
 10   format(30i4)
      write(5,10)ist,ist1,imin,imin1
      write(5,*)best,best1
      do 11 i=1,nn
      write(5,10)(itab(i,j),j=1,kk),(ifr(i,j),j=1,kk)
     * ,(ifrs(i,j),j=1,kk)
 11   continue
      open(6,file='inp')
      do 50 j=1,lines
      read(6,1)(lst(j,k),k=1,60)
 50   continue
      close(6)
      open(7,file='prog.lis')
      do 7 i=1,lines
      write(7,1)(lst(i,k),k=1,60)
 7    continue
      do 8 i=1,nn
      do 8 j=1,kk
      if(ifr(i,j).eq.0)go to 8
      write(7,9)i,j
 9    format('FR LX(',i2,',',i2,')')
 8    continue
      write(7,*)'OU ul ND=5'
      close(7)
      stop
 2    read(5,10)ist,ist1,imin,imin1
      read(5,*)best,best1
      do 12 i=1,nn
      read(5,10)(itab(i,j),j=1,kk),(ifr(i,j),j=1,kk)
     * ,(ifrs(i,j),j=1,kk)
 12   continue
      close(5)
      open(8,file='prog.out')
 21   read(8,20)(lin(i),i=1,80)
 20   format(80a1)
      do 60 ik=5,15
      if(lin(ik+1).ne.numb(12))go to 60
      if(lin(ik+2).ne.numb(13))go to 60
      if(lin(ik+3).ne.numb(14))go to 60
      go to 61
 60   continue
      go to 21
 61   continue
      idf=find(lin,60,numb)
      chi=find(lin(30),38,numb)
      ff=chi/idf
      write(0,*)ff,best,best1
      if(ff.lt.best)go to 30
      if(ff.lt.best1)go to 31
      go to 35
 30   best=ff
      do 40 i=1,nn
      do 40 j=1,kk
 40   itab(i,j)=0
      open(9,file='best.lis')
      open(7,file='prog.lis')
 33   read(7,1,end=32)(lin(i),i=1,80)
      write(9,1)(lin(i),i=1,80)
      go to 33
 32   close(7)
      close(9)
 38   ist=1
      ist1=1
      best1=100000.
 36   ifr(ist,ist1)=1-ifr(ist,ist1)
      go to 99
 31   if(itab(ist,ist1).gt.0.and.iter-itab(ist,ist1).lt.itmax*nn*kk)
     * go to 35
      do 34 i=1,nn
      do 34 j=1,kk
 34   ifrs(i,j)=ifr(i,j)
      imin=ist
      imin1=ist1
      best1=ff
 35   ifr(ist,ist1)=1-ifr(ist,ist1)
      ist=ist+1
      if(ist.le.nn)go to 36
      ist=1
      ist1=ist1+1
      if(ist1.le.kk)go to 36
      do 37 i=1,nn
      do 37 j=1,kk
 37   ifr(i,j)=ifrs(i,j)
      itab(imin,imin1)=iter
      go to 38
      end
      function find(lin,n,numb)
      character*1 lin(60),numb(14)
      find=0.
      idiv=0
      fact=1.
      fact1=10.
      do 1 i=1,n
      do 2 j=1,11
      jj=j
      if(lin(i).eq.numb(j))go to 3
  2   continue
      if(find.eq.0.)go to 1
      return
  3   if(jj.le.10)go to 4
      idiv=1
      fact=0.1
      fact1=1.
      go to 1
  4   find=find*fact1+fact*(jj-1)
      if(idiv.eq.1)fact=fact*0.1
  1   continue
      end
