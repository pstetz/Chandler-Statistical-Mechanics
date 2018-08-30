10 REM PROGRAM TO SIMULATE A TWO-LEVEL SYSTEM  COUPLED TO
20 REM AN ADIABATIC GAUSSIAN FLUCTUATING FIELD
30 REM version#1 on DEC 09 85 by FARAMARZ RABII at U of Penn
40 ON KEY(2) GOSUB 80
50 KEY(2) ON
60 DIM SIGMA%(32)
70 DIM ESOL!(2000)
80 REM RESTART FROM HERE
90 CLS
100 GOSUB 1330
110 PRINT"Would you like an introduction to the methods used here ?"
120 INPUT"*** YES, or NO ***";K$
130 IF (K$<>"yes") AND (K$<>"YES") AND (K$<>"no") AND (K$<>"NO") THEN GOTO 110
140 IF (K$="yes") OR (K$="YES") THEN GOSUB 2190
150 CLS
160 RANDOMIZE(310941!)
170 KEY(2) ON
180 LOCATE 1,1
190 PRINT "ENTER THE LOCALIZATION PARAMETER 'L'"
200 INPUT "***BETWEEN 0.01 AND 10 ***";LOCAL!
210 IF LOCAL!<.01 GOTO 200
220 IF LOCAL!>10 GOTO 200
230 PRINT "ENTER THE REDUCED TEMPERATURE 'BETA'";
240 INPUT "***BETWEEN 0.01 AND 16 ***";RBETA!
250 IF RBETA!<.01 GOTO 240
260 IF RBETA!>16 GOTO 240
270 INPUT "TOTAL NUMBER OF STEPS (INTEGER BELOW 1E+20)";MOVE!
280 IF MOVE!<0 THEN GOTO 270
290 IF MOVE!<>INT(MOVE!) THEN GOTO 270
300 IF MOVE!>1E+20 THEN GOTO 270
310 IVAR!=RBETA!/(2*LOCAL!)
320 FSTEP!=IVAR!
330 KAY!=-.5*LOG(RBETA!/32)
340 TKAY!=KAY!*2
350 LAMBDA!=RBETA!/32
360 TLAMBDA!=2*LAMBDA!
370 NEWENERGY!=0
380 FILD!=0
390 ESOLV!=0
400 AVDIP!=0
410 FOR I=1 TO 32
420 IF I/2=INT(I/2) THEN SIGMA%(I)=1 ELSE SIGMA%(I)=-1
430 NEXT I
440 CLS
450 GOSUB 1680
460 GOSUB 1070
470 GOSUB 950
480 FOR IRUN!=1 TO MOVE!
490 LOCATE 16,45
500 PRINT "STEP =";IRUN!
510 JAY%=INT(RND*39)+1
520 IF JAY%>32 THEN GOTO 560
530 OLDSIGMA%=SIGMA%(JAY%)
540 IF OLDSIGMA%=1 THEN SIGMA%(JAY%)=-1 ELSE SIGMA%(JAY%)=1
550 GOTO 590
560 GOSUB 1560
570 OLDFILD!=FILD!
580 FILD!=FILD!+FSTEP!*(.5-RND)
590 OLDENERGY!=NEWENERGY!
600 OLDESOLV!=NEWESOLV!
610 GOSUB 1140
620 IF NEWENERGY!<OLDENERGY! THEN GOTO 670
630 REM accept move unit probability
640 IF JAY%<33 THEN DIP!=NEWDIP!
650 GOSUB 1450
660 GOTO 760
670 PROB!=EXP(NEWENERGY!-OLDENERGY!)
680 IF RND>PROB! THEN GOTO 720
690 IF JAY%<33 THEN DIP!=NEWDIP!
700 GOSUB 1450
710 GOTO 760
720 IF JAY%>32 THEN FILD!=OLDFILD! ELSE SIGMA%(JAY%)=OLDSIGMA%
730 NEWENERGY!=OLDENERGY!
740 NEWESOLV!=OLDESOLV!
750 GOSUB 1450
760 ESOLV!=ESOLV!+NEWESOLV!
770 TST!=IRUN!/100
780 IF TST!<>INT(TST!) THEN GOTO 920
790 ESOL!(TST!)=ESOLV!/IRUN!
800 LOCATE 17,22
810 PRINT "               "
820 LOCATE 17,3
830 PRINT "SOLVATION ENERGY = ";-NEWESOLV!/RBETA!
840 LOCATE 18,3
850 PRINT "(IN UNITS OF DELTA)"
860 LOCATE 19,14
870 PRINT "               "
880 LOCATE 19,3
890 PRINT"DIPOLE = ";DIP!/32
900 LOCATE 20,3
910 PRINT "(IN UNITS OF MU)"
920 IF IRUN!>MOVE! THEN GOTO 2010
930 NEXT IRUN!
940 GOTO 2010
950 REM OUT PUT CONFIGURATION
960 DIP!=0
970 LOCATE 11,3
980 FOR INDEX%=1 TO 32
990 IF SIGMA%(INDEX%)=1 THEN PRINT CHR$(176);CHR$(176); ELSE PRINT CHR$(219);CHR$(219);
1000 DIP!=DIP!+SIGMA%(INDEX%)
1010 NEXT INDEX%
1020 LOCATE 12,3
1030 FOR INDEX%=1 TO 32
1040 IF SIGMA%(INDEX%)=1 THEN PRINT CHR$(219);CHR$(219); ELSE PRINT CHR$(176);CHR$(176);
1050 NEXT INDEX%
1060 RETURN
1070 REM subsection to compute overall energy
1080 NEWENERGY!=0
1090 FOR I=1 TO 32
1100 IF I=32 THEN J=1 ELSE J=I+1
1110 NEWENERGY!=NEWENERGY!+KAY!*SIGMA%(I)*SIGMA%(J)
1120 NEXT I
1130 RETURN
1140 REM SUBSECTION TO QUICKLY COMPUTE NEW ENERGY
1150 IF JAY%=1 THEN JAYM1%=32 ELSE JAYM1%=JAY%-1
1160 IF JAY%=32 THEN JAYP1%=1 ELSE JAYP1%=JAY%+1
1170 IF JAY%>32 THEN GOTO 1290
1180 IF SIGMA%(JAY%)=1 THEN GOTO 1240
1190 NEWESOLV!=OLDESOLV!-FILD!*TLAMBDA!
1200 NEWENERGY!=OLDENERGY!-OLDESOLV!+NEWESOLV!
1210 NEWENERGY!=NEWENERGY!-TKAY!*(SIGMA%(JAYP1%)+SIGMA%(JAYM1%))
1220 NEWDIP!=DIP!-2
1230 GOTO 1320
1240 NEWESOLV!=OLDESOLV!+FILD!*TLAMBDA!
1250 NEWENERGY!=OLDENERGY!-OLDESOLV!+NEWESOLV!
1260 NEWENERGY!=NEWENERGY!+TKAY!*(SIGMA%(JAYP1%)+SIGMA%(JAYM1%))
1270 NEWDIP!=DIP!+2
1280 GOTO 1320
1290 NEWENERGY!=OLDENERGY!+IVAR!*(OLDFILD!*OLDFILD!-FILD!*FILD!)
1300 NEWESOLV!=LAMBDA!*DIP!*FILD!
1310 NEWENERGY!=NEWENERGY!-OLDESOLV!+NEWESOLV!
1320 RETURN
1330 REM SUBSECTION FOR INTRODUCTION
1340 PRINT "          SIMULATION OF A TWO-LEVEL SYSTEM COUPLED"
1350 PRINT "                    TO A GAUSSIAN FIELD"
1360 PRINT
1370 PRINT "                            BY:"
1380 PRINT
1390 PRINT "                       FARAMARZ RABII"
1400 PRINT
1410 PRINT"     This  program simulates  a quantum mechanical one dimensional"
1420 PRINT"dipole  coupled to  an  adiabatically  fluctuating  field  obeying"
1430 PRINT"gaussian  statistics."
1440 RETURN
1450 REM update display
1460 IF (IRUN!/4)<>INT((IRUN!/4)) THEN RETURN
1470 IF JAY%>32 THEN GOTO 1520
1480 LOCATE 11,2*JAY%+1
1490 IF SIGMA%(JAY%)=1 THEN PRINT CHR$(176);CHR$(176) ELSE PRINT CHR$(219);CHR$(219)
1500 LOCATE 12,2*JAY%+1
1510 IF SIGMA%(JAY%)=1 THEN PRINT CHR$(219);CHR$(219) ELSE PRINT CHR$(176);CHR$(176)
1520 RETURN
1530 LOCATE INDEX%+12,75
1540 PRINT " "
1550 NEXT INDX%
1560 REM SUBSECTION TO OUT-PUT EXTERNAL FIELD
1570 OLDFLDMAG%=FLDMAG%
1580 FLDMAG%=-INT(FILD!/FSTEP!)
1590 IF ABS(FLDMAG%)>9 THEN FLDMAG%=9*SGN(FLDMAG%)
1600 IF ABS(FLDMAG%)>=ABS(OLDFLDMAG%) THEN GOTO 1640
1610 LOCATE OLDFLDMAG%+11,74
1620 PRINT CHR$(219);CHR$(219);CHR$(219)
1630 GOTO 1670
1640 LOCATE FLDMAG%+11,74
1650 PRINT CHR$(176);CHR$(176);CHR$(176)
1660 GOTO 1670
1670 RETURN
1680 REM SUBSECTION TO DRAW A BOX AROUND THE FIELD OUT-PUT
1690 FOR I=1 TO 19
1700 LOCATE I+1,73
1710 PRINT CHR$(179);CHR$(219);CHR$(219);CHR$(219);CHR$(179)
1720 NEXT I
1730 LOCATE 1,73
1740 PRINT "-----"
1750 LOCATE 21,73
1760 PRINT "-----"
1770 LOCATE 11,70
1780 PRINT "0.0"
1790 LOCATE 2,70
1800 PRINT "9.0"
1810 LOCATE 20,69
1820 PRINT "-9.0"
1830 LOCATE 11,74
1840 PRINT CHR$(176);CHR$(176);CHR$(176)
1850 LOCATE 22,75
1860 PRINT CHR$(24)
1870 LOCATE 23,63
1880 PRINT "FLUCTUATING FIELD";
1890 LOCATE 23,9
1900 PRINT "TO RESTART PRESS F2"
1910 REM DISPLAY INITIAL CONDITIONS
1920 LOCATE 3,1
1930 PRINT "LOCALIZATION PARAMETER L=";LOCAL!
1940 PRINT
1950 PRINT "REDUCED TEMPERATURE BETA=";RBETA!
1960 LOCATE 9,29
1970 PRINT "THE QUANTUM PATH"
1980 LOCATE 10,36
1990 PRINT CHR$(25)
2000 RETURN
2010 REM OUT-PUT SOLVATION ENERGY VALUES
2020 LPRINT "RESULTS FOR THE SIMULATION OF A TWO-LEVEL SYSTEM COUPLED"
2030 LPRINT "TO AN ADIABATIC FIELD."
2040 LPRINT
2050 LPRINT "INITIAL CONDITIONS ARE:"
2060 LPRINT
2070 LPRINT "LOCALIZATION PARAMETER  L = ";LOCAL!
2080 PRINT
2090 LPRINT "REDUCED TEMPERATURE BETA = ";RBETA!
2100 LPRINT
2110 LPRINT "#OF STEPS","AVERAGE SOLVATION ENERGY IN UNITS OF DELTA"
2120 LPRINT
2130 IMAX%=MOVE!/100
2140 FOR I=1 TO IMAX%
2150 LPRINT I*100,-ESOL!(I)/RBETA!
2160 NEXT I
2170 LOCATE 23,1
2180 END
2190 CLS
2200 PRINT
2210 PRINT"                          INTRODUCTION"
2220 PRINT
2230 PRINT"     The  goal  of  this program is to compute equilibrium averages"
2240 PRINT"for  a  system  composed  of a quantum-mechanical dipole coupled to"
2250 PRINT"an  adiabatically  fluctuating  field . The  dipole is modeled as a"
2260 PRINT"two level  system  with  the  states  corresponding  to  the dipole"
2270 PRINT"pointing  along  or  opposite  the  field . The  partition function"
2280 PRINT"is  then  written  as  the  trace of  the  density matrix using the"
2290 PRINT"Feynman  path integrals . The  trace  becomes  a path integral over"
2300 PRINT"the  field  and  a path sum over the two level system's states."
2310 PRINT"     Since  the  field is adiabatic, the path integral reduces to a"
2320 PRINT"single  integration . The  two level system's path sum (technically"
2330 PRINT"an  infinite  number  of summations ) is  approximated  as a finite"
2340 PRINT"number  of  summations ( 32 here ) . The  problem is now isomorphic"
2350 PRINT"to  that  of  a  one  dimensional  classical ring  of  32  electric"
2360 PRINT"dipoles  with  nearest   neighbor  interactions  in  a  fluctuating"
2370 PRINT"external  electric  field . Each term in the summations corresponds"
2380 PRINT"to a configuration of the ring which we call 'The Quantum Path'."
2390 PRINT
2400 INPUT"PRESS RETURN WHEN READY";K$
2410 CLS
2420 PRINT
2430 PRINT"                PARAMETERIZATION OF THE SYSTEM"
2440 PRINT
2450 PRINT"     The  system adressed here depends on the following parameters:"
2460 PRINT"1)Mu, the  magnitude of the dipole ( diagonal matrix element of the"
2470 PRINT"two level system)."
2480 PRINT"2)Delta, half  the  tunnel  splitting  of the dipole ( off diagonal"
2490 PRINT"matrix element of the two level system )."
2500 PRINT"3)Sigma, the variance of the gaussian field."
2510 PRINT"4)T, the temperature."
2520 PRINT"     Although , there are four parameters above, there are only two"
2530 PRINT"energies   in  the  problem : one  corresponding  to  the tunneling"
2540 PRINT"frequency  of  the  dipole , and  the  other  corresponding  to the"
2550 PRINT"interaction with the field. This makes it possible to formulate the"
2560 PRINT"problem in terms of two unitless reduced parameters given below:"
2570 PRINT"1)The localization parameter L = (sigma*mu*mu/delta)"
2580 PRINT"2)The reduced temprature beta = delta/(Boltzmann constant*T) "
2590 PRINT
2600 INPUT"PRESS RETURN WHEN READY";K$
2610 CLS
2620 PRINT"                      RUNNING THE PROGRAM"
2630 PRINT
2640 PRINT"     When L < 1 , adiabatic  electric potentials in both ground and"
2650 PRINT"excited  states  are  monostable . But when L > 1, the ground state"
2660 PRINT"potential  as  a  function  of  the  field  exhibits two degenerate"
2670 PRINT"minima . This  behavior  leads to quenching of the dipole resonance"
2680 PRINT"or tunneling. The effect becomes a singular one at zero temperature"
2690 PRINT"( beta  goes  to  infinity ) where there  is a symmetry breaking of"
2700 PRINT"the quantal dipole."
2710 PRINT"     The  program  will ask you to enter the parameters L and beta."
2720 PRINT"Once  these  parameters  are entered, the program will simulate the"
2730 PRINT"system , computing  the  average  interaction  energy  of the field"
2740 PRINT"with  the  dipole ( in  units  of delta ) and  storing it every 100"
2750 PRINT"steps. After the  prescribed  number  of  steps , the  interaction"
2760 PRINT"energy recorded  100 steps is printed out."
2770 PRINT"     The  initial  state  of the  two-level  system  has the dipole"
2780 PRINT"resonating between it's two configurations. The display illustrates"
2790 PRINT"the  quantum path for the dipole showing the progress of the Monte-"
2800 PRINT"Carlo trajectory  for  the  32 sections of the quantum path and the"
2810 PRINT"the fluctuating field in units of delta/mu . Every  100  steps  the"
2820 PRINT"display  updates  the  instantaneous  value  of   the  field/dipole"
2830 PRINT"interaction energy. Updates of the display are made every 4 steps."
2840 PRINT
2850 INPUT "PRESS RETURN WHEN READY";K$
2860 RETURN
2900 END
