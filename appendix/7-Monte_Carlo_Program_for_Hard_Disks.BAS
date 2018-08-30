REM***MONTE CARLO
REM***an input file "data" is needed where the first n lines contain
REM***the coordinates of the disks, the next line contains the total
REM***number of passes, and the last mg lines the total number found
REM***in each g(r) shell. "data" is updated with each run.
REM***"del" is the maximum step size.  "n" is the number of disks.
REM***"rho" is the density of disks.  The screen is printed every "mcopy"
REM***times.  "mg" is the number of divisions in G(r).  "pas" is the number
REM*** of passes.  "rgmax" is the maximum of the range G(r) is calculated
REM*** over.  G(r=1) is found by linear extrapolation.  ---John McCoy---
REM*** To stop program, type "command ." followed by "menu reset"
del=.05
n=20
mcopy=1000
mg=15
PRINT "Input number of passes."
INPUT pas
rgmax=2.5
dgr=(rgmax-1)/mg
pac = 100/n
PRINT " If you wish to start from previous run in file data, type 1;"
PRINT " If you wish to start new run, type 0."
INPUT istart
IF istart=0 THEN PRINT "Input a density between 0 and 0.85."
IF istart=0 THEN INPUT rho
IF istart=1 THEN PRINT "input 1 if you wish to keep old g(r) data, else input 0."
igr=1
IF istart=1 THEN INPUT igr
DIM f(15),fmn(15),r(15),ddr(15),g(15)
IF istart=1 THEN OPEN "data" FOR INPUT AS #1
IF istart=0 THEN OPEN "datazero" FOR INPUT AS #1
IF istart=1 THEN INPUT #1, rho
ppres =rho*3.141593/2
sig=(rho/n)^.5
rmax=4*sig
srgmax=rgmax*sig
rcir%=FIX(100*sig)

FOR i=1 TO mg
r(i)=1+(i-.5)*dgr
NEXT i

dd=3.14*(n-1)*rho
FOR i=1 TO mg
ddr(i)=dd*((1+(i)*dgr)^2-(1+(i-1)*dgr)^2)
NEXT i

t$="monte carlo"
IN:
CALL HIDECURSOR
LET today$=DATE$
MENU 1,0,1," Monte Carlo            "
MENU 2,0,1," "
MENU 3,0,1,"Hard Disks "
MENU 4,0,1,"       "
MENU 5,0,1,today$
WINDOW 1,t$,(0,20)-(550,350),3
DIM x(20),y(20)
RANDOMIZE TIMER

FOR I=1 TO n
INPUT #1,x(I),y(I)
NEXT I
INPUT #1,pps
IF igr=0 THEN pps=0
FOR i=1 TO mg
INPUT #1,f(i)
IF igr=0 THEN f(i)=0
NEXT i
CLOSE #1

drw:
CALL MOVETO(240,25)
PRINT 5,
CALL MOVETO(240,65)
PRINT 4,
CALL MOVETO(240,105)
PRINT 3,
CALL MOVETO(240,145)
PRINT 2,
CALL MOVETO(240,185)
PRINT 1,
CALL MOVETO(240,225)
PRINT 0,
CALL MOVETO(240,125)
G$="g"
PRINT G$,
CALL MOVETO(259,240)
PRINT 1,
CALL MOVETO(339,240)
PRINT 2,
CALL MOVETO(419,240)
PRINT 3,
CALL MOVETO(359,255)
r$="r"
PRINT r$,

CALL MOVETO (20,20)
CALL LINE (0,200)
CALL LINE (200,0)
CALL LINE (0,-200)
CALL LINE (-200,0)
CALL MOVETO(270,20)
CALL LINE (0,200)
CALL LINE(160,0)
CALL MOVETO(270,20)
CALL LINE (-5,0)
CALL MOVETO(270,60)
CALL LINE (-5,0)
CALL MOVETO(270,100)
CALL LINE (-5,0)
CALL MOVETO(270,140)
CALL LINE (-5,0)
CALL MOVETO(270,180)
CALL LINE (-5,0)
CALL MOVETO(270,220)
CALL LINE (-5,0)
CALL MOVETO(270,220)
CALL LINE (0,5)
CALL MOVETO(310,220)
CALL LINE (0,5)
CALL MOVETO(350,220)
CALL LINE (0,5)
CALL MOVETO(390,220)
CALL LINE (0,5)
CALL MOVETO(430,220)
CALL LINE (0,5)
CALL MOVETO(40,250)
pas$="passes="
PRINT pas$,pps,
CALL MOVETO(40,265)
mov$="particles="
PRINT mov$,n,
CALL MOVETO(40,280)
den$="density="
PRINT den$,rho,
g1$ = "g(r=1) ="
pre$ ="P§/rho ="
ac$ = "%accept="


FOR i=1 TO n
xx=200*x(i)+20
yy=200*y(i)+20
CIRCLE (xx,yy),rcir%,33
NEXT i

FOR k=1 TO pas

acc% = 0

FOR j=1 TO n

r=1-2*RND(1)
xn=x(j)+del*r
r=1-2*RND(1)
yn=y(j)+del*r
xxo1=x(j)
yyo1=y(j)
xxo=200*xxo1+20
yyo=200*yyo1+20

FOR jj=1 TO mg
fmn(jj)=f(jj)
NEXT jj

FOR ij=1 TO n
IF ij=j THEN GOTO 10
rx= x(ij)-xn
ry= y(ij)-yn
IF rx>.5 THEN rx=rx-1 ELSE IF rx<-.5 THEN rx=rx+1
IF ry>.5 THEN ry=ry-1 ELSE IF ry<-.5 THEN ry=ry+1
r=(rx^2+ry^2 )^.5
IF r<sig THEN GOTO new1
IF r>srgmax  THEN GOTO 10
xxx=((r/sig)-1)/dgr
ii=FIX(xxx)+1
fmn(ii)=fmn(ii)+1 


10 :  
NEXT ij

acc% = acc% + 1

GOTO new2

 NEW1:
xn =xxo1
yn=yyo1

FOR jj=1 TO mg
fmn(jj)=f(jj)
NEXT jj

FOR ij=1 TO n
IF ij=j THEN GOTO 20
rx= x(ij)-xn
ry= y(ij)-yn
IF rx>.5 THEN rx=rx-1 ELSE IF rx<-.5 THEN rx=rx+1
IF ry>.5 THEN ry=ry-1 ELSE IF ry<-.5 THEN ry=ry+1
r=(rx^2+ry^2 )^.5
IF r>srgmax  THEN GOTO 20
xxx=((r/sig)-1)/dgr
ii=FIX(xxx)+1
fmn(ii)=fmn(ii)+1 

20 :
NEXT ij

NEW2:
FOR jj=1 TO mg
f(jj)=fmn(jj)
NEXT jj

x(j)=xn
y(j)=yn
IF x(j)<0 THEN x(j)=x(j)+1 ELSE IF x(j)>1 THEN x(j)=x(j)-1
IF y(j)<0 THEN y(j)=y(j)+1 ELSE IF y(j)>1 THEN y(j)=y(j)-1
xx=200*x(j)+20
yy=200*y(j)+20
CIRCLE (xxo,yyo),rcir%,30
CIRCLE (xx,yy),rcir%,33

NEXT j

CALL MOVETO (20,20)
CALL LINE (0,200)
CALL LINE (200,0)
CALL LINE (0,-200)
CALL LINE (-200,0)
pps=pps+1

FOR jj=1 TO mg
g(jj)=f(jj)/(ddr(jj)*pps)
NEXT jj

CALL MOVETO((270+(r(1)-1)*80),(220-g(1)*40))
rec%(0)=0
rec%(1)=271
rec%(2)=219
rec%(3)=550
CALL ERASERECT (VARPTR(rec%(0)))

FOR jj=1 TO mg-1
xg=(r(jj+1)-r(jj))*80
yg=(g(jj)-g(jj+1))*40
CALL LINE(xg,yg)
NEXT jj
gcont = 1.5*g(1) - .5*g(2)
pres = 1 + ppres*gcont
CALL MOVETO(240,280)
PRINT g1$,gcont,
CALL MOVETO(240,295)
PRINT pre$,pres,
pacc = acc%*pac
CALL MOVETO(40,295)
PRINT ac$,pacc,

CALL MOVETO(40,250)
PRINT pas$,pps,
pptest=mcopy*FIX(pps/mcopy)
IF pptest=FIX(pps) THEN LCOPY

NEXT k

 out:
 
OPEN "data" FOR OUTPUT AS #2
WRITE #2,rho
FOR I =1 TO n
xn=10000*x(i)
yn=10000*y(i)
xn=FIX(xn)/10000
yn=FIX(yn)/10000
WRITE #2,xn,yn
NEXT i
 
WRITE #2,pps
 
FOR i=1 TO mg
WRITE #2,f(i)
NEXT i
CLOSE #2
WINDOW CLOSE 1
MENU RESET
END
