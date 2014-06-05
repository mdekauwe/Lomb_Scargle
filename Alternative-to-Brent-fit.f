
!   Test for best AR(1) spectrum fit to the log of the median
!   smoothed spectrum by varying Rho and SO (via least squares).
!   Regular-spaced data: use lowest of either Rho (initial estimate)
!   or Rho (Mudelsee method) (=ROH) as upper bound for rho during fitting
!   (i.e. vary ROH and SOFA).
!   Irregularly-spaced data: use Rho (Mudelsee method) (i.e. vary
!   SOFA only).
       IF(NFIT == 0.OR.NFIT == 1)THEN
!   If on Monte Carlo spectral run miss out screen notes.
         IF(NSPECB == 0)THEN
           WRITE(*,*)' ******************************************'
           WRITE(*,*)' *** Spectral background location via:  ***'
           IF(MUDEST == 0)THEN
             IF(NFIT == 0)THEN
               WRITE(*,*)' *** White Noise model conf. levels.    ***'
             ENDIF
             IF(NFIT == 1)THEN
               WRITE(*,*)' *** a robust-fitted AR(1) model        ***'
               WRITE(*,*)' *** for finding confidence levels      ***'
               WRITE(*,*)' *** (Mann & Lees, 1996, Clim. Change). ***'
             ENDIF
           ELSE
             WRITE(*,*)' *** fitting of AR(1) model based on    ***'
             WRITE(*,*)' *** best estimate of the lag-1         ***'
             WRITE(*,*)' *** correlation coefficient & finding  ***'
             WRITE(*,*)' *** best-fit mean log spectrum level   ***'
             WRITE(*,*)' *** in to find confidence levels.      ***'
           ENDIF
           WRITE(*,*)' ******************************************'
         ENDIF
         IF((XLOUT(NMED)+XLOUT(2)) == 0.0D0)THEN
           WRITE(*,*)
           WRITE(*,*)' *** Problem initialising mid spectrum ***'
           STOP
         ENDIF
!   Find maximum and minimum log median-smoothed power value (in XLOUT)
         XMAXLP=XLOUT(2)
         XMINLP=XLOUT(NMED) 
         DO 335 I=2,NMED
           IF(XLOUT(I) > XMAXLP)XMAXLP=XLOUT(I)
           IF(XLOUT(I) < XMINLP)XMINLP=XLOUT(I)
  335    ENDDO
         SOLMID=(XMAXLP+XMINLP)/2.0D0
         RANGE=(XMAXLP-XMINLP)/2.0D0
         IF(NFIT == 0)RANGE=XMAXLP-SOLMID
         SOLMAX=SOLMID+RANGE
         SOLMIN=SOLMID-RANGE
         SOF=0.0D0
         SOLINT=0.01D0
         IF(ROH <= 0.0D0)NFIT=0
!   Using AR1 fit and  uniformly/regularly-spaced data vary ROH and SOFA.
         IF(NFIT == 1.AND.NIRR == 0)THEN
           SOFIT=SOLMAX+SOLINT
           ROHB=ROH
           XLSIN=1000000000000.0
           XLSMIN=XLSIN
  338      FORMAT('           Testing LogSO = ',F7.4,' to ',F7.4)
  339      FORMAT('                     ROH =  0.000  to  ',F5.3)
           WRITE(*,*)
           WRITE(*,338)SOLMAX,SOLMIN
           WRITE(*,339)ROH
           WRITE(*,*)
!   Try repeatedly reducing SOFIT to find a better fit:
           CTLIM=0
           DO 340
             CTLIM=CTLIM+1
             IF(CTLIM > N)THEN
               WRITE(*,*)' *** EXIT DO 340 CTLIM>N ***'
             ENDIF
             IF(CTLIM > N)EXIT                 ! EXIT DO 340
             SOFIT=SOFIT-SOLINT    
             IF(SOFIT < SOLMIN)EXIT            ! EXIT DO 340
             SOFA=10.0**SOFIT
!   Call minimization procedure (using functions BRENT2 and XLS2) to find
!   best-fit value for rho given current SOFA.
             GUESS=ROHB-0.1D0
             IF(GUESS < 0.0D0)GUESS=ROHB/2.0D0
! Old linear freq:  DUM4=BRENT2(0.0D0,GUESS,ROHB,XLS2,TOL,ROHT,XLOUT,XFR,
!     *       NMED,SOFA,XNYQ)
             XLNYQ=LOG10(XNYQ)
             DUM4=BRENT2(0.0D0,GUESS,ROHB,XLS2,TOL,ROHT,XLOUT,XLFR,
     *       NMED,SOFA,XLNYQ)
! Old:              DUM4=BRENT2(0.0D0,GUESS,ROHB,XLS2,TOL,ROHT,XOUT,XFR,
!     *       NMED,SOFA,XNYQ)
!   Test for minimum least squares given current value of SOFA
             IF(ROHT /= 2.0D0)THEN
               IF(ROHT < 0.0D0)THEN
                 WRITE(*,*)' ROHT <= 0.0'
                 EXIT                          ! EXIT DO 340
               ENDIF
               IF(ROHT > ROHB)THEN
                 WRITE(*,*)' ROHT > ROHB'
                 EXIT                          ! EXIT DO 340
               ENDIF
               IF(DUM4 < XLSMIN)THEN
                 XLSMIN=DUM4
                 ROHF=ROHT
                 SOF=SOFA
                 CYCLE                         ! CYCLE DO 340
               ENDIF
               IF(DUM4 /= XLSIN)CYCLE          ! CYCLE DO 340
             ENDIF
             EXIT                              ! EXIT DO 340
  340      ENDDO
           IF(XLSMIN == XLSIN)THEN
             WRITE(*,*)' *** Warning: Least Sqs min = Init. value! ***'
             STOP
           ENDIF
           IF(ROHT < 0.0D0.OR.ROHT > 1.0D0)THEN
             WRITE(*,*)' *** Problem with Rho fitting ***'
             STOP
           ENDIF
         ENDIF
!   Find best fit value of SOF given ROH (for regularly-spaced data
!   and white noise fit or irregularly spaced data).
         IF(NFIT == 0.OR.(NIRR == 1.AND.NFIT == 1))THEN
           SOLMI=10.0**SOLMIN
           SOLMA=10.0**SOLMAX
           SOLMD=10.0**SOLMID
  350      FORMAT('           Testing LogSO= ',F7.4,' to ',F7.4)
  351      FORMAT('               using ROH=  ',F5.3)
           IF(ROH <= 0.0D0)ROH=0.0D0
           ROHB=ROH
           IF(NFIT == 0)ROHB=0.0D0
           WRITE(*,350)SOLMAX,SOLMIN
           WRITE(*,351)ROHB
           WRITE(*,*)
!   Call minimization procedure (using functions BRENT3 and XLS3) to find
!   best-fit value for SOFA given ROH.
! Old:       DUM5=BRENT3(SOLMI,SOLMD,SOLMA,XLS3,TOL,SOFT,XLOUT,
!     *     XFR,NMED,ROHB,XNYQ)
           XLNYQ=LOG10(XNYQ)
           DUM5=BRENT3(SOLMI,SOLMD,SOLMA,XLS3,TOL,SOFT,XLOUT,
     *     XLFR,NMED,ROHB,XLNYQ)
! Old:       DUM5=BRENT3(SOLMI,SOLMD,SOLMA,XLS3,TOL,SOFT,XOUT,
!     *     XFR,NMED,ROHB,XNYQ)
!   Output is least squares best fit of SOF.
           ROHF=ROHB
           IF(SOFT == 20.0D0)THEN
             SOF=SOLMID
           ELSE
             SOF=SOFT
           ENDIF
         ENDIF
!   Check BRENT result for false XLS minimum.
         NPOS=0
         NNEG=0
         XLSUM=0.0d0
         DO 362 I=1,NMED
           XARP2=(SOF*(1.0-(ROHF*ROHF)))/(1.0-2*ROHF*COS((XFR(I)/XNYQ)*
     *     XPI)+(ROHF*ROHF))
           XLP2=LOG10(XARP2)
           XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))
           IF(XLP2-XLOUT(I) > 0.0)NPOS=NPOS+1
           IF(XLP2-XLOUT(I) < 0.0)NNEG=NNEG+1
  362    ENDDO
         BRENTSUM=XLSUM
!   If BRENT fitted background spectrum (almost) entirely above median smoothed
!   spectrum then drop fitted background to find true XLS minimum.
         XLMIN=BRENTSUM 
         IF(SOF <= 0.0)THEN
           WRITE(*,*)' *** SOF negative or zero! ***'
           STOP
         ENDIF      
         XLSOF=LOG10(SOF)
         SOF2=0.0
         PROP=NPOS/NMED
         IF(PROP > 0.7)THEN
           DO 364
             XLSOF=XLSOF-0.1
             IF(XLSOF < SOLMIN)EXIT          ! EXIT DO 364
             SOF2=10.0**XLSOF
             XLSUM=0.0d0
             DO 363 I=1,NMED
               XARP2=(SOF2*(1.0-(ROHF*ROHF)))/(1.0-2*ROHF*COS((XFR(I)/
     *         XNYQ)*XPI)+(ROHF*ROHF))
               XLP2=LOG10(XARP2)
               XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))             
  363        ENDDO
             IF(XLSUM < XLMIN)THEN
               XLMIN=XLSUM
               SOFNEW=SOF2
             ENDIF
             CYCLE                           ! CYCLE DO 364
  364      ENDDO  
         ENDIF
!   If BRENT fitted background spectrum (almost) entirely below median smoothed
!   spectrum then raise fitted background to find true XLS minimum.
         PROP=NNEG/NMED
         IF(PROP > 0.7)THEN
           DO 366
             XLSOF=XLSOF+0.1
             IF(XLSOF > SOLMAX)EXIT          ! EXIT DO 366
             SOF2=10.0**XLSOF
             XLSUM=0.0d0
             DO 365 I=1,NMED
               XARP2=(SOF2*(1.0-(ROHF*ROHF)))/(1.0-2*ROHF*COS((XFR(I)/
     *         XNYQ)*XPI)+(ROHF*ROHF))
               XLP2=LOG10(XARP2)
               XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))             
  365        ENDDO
             IF(XLSUM < XLMIN)THEN
               XLMIN=XLSUM
               SOFNEW=SOF2
             ENDIF
             CYCLE                           ! CYCLE DO 366
  366      ENDDO  
         ENDIF
         IF(XLMIN < BRENTSUM)THEN
           SOF=SOFNEW
           WRITE(*,*)' Spectral background from fixing BRENT fit'
         ENDIF
!   If final fitted ROH (= ROHF) is not raw or Mudelesee ROH then 
!   check simple alternative to fit from BRENT output by using fixed
!   initial ROH (=ROHI).
         XLREV=XLMIN
         IF(ROHF /= ROHI)THEN
           XLSOF=SOLMAX+0.1
           DO 368 
             XLSOF=XLSOF-0.1
             IF(XLSOF < SOLMIN)EXIT          ! EXIT DO 368
             SOF2=10.0**XLSOF
             XLSUM=0.0d0
             DO 367 I=1,NMED
               XARP2=(SOF2*(1.0-(ROHI*ROHI)))/(1.0-2*ROHI*COS((XFR(I)/
     *         XNYQ)*XPI)+(ROHI*ROHI))
               XLP2=LOG10(XARP2)
               XLSUM=XLSUM+((XLP2-XLOUT(I))*(XLP2-XLOUT(I)))             
  367        ENDDO
             IF(XLSUM < XLMIN)THEN
               XLMIN=XLSUM
               SOFNEW=SOF2
             ENDIF
             CYCLE                          ! CYCLE DO 368             
  368      ENDDO
         ENDIF
         IF(XLMIN < XLREV)THEN
           SOF=SOFNEW
           ROHF=ROHI
           WRITE(*,*)' Spectral background from ROH raw gives best fit'
         ENDIF
         IF(SOF > 0.0D0)XLSOF=LOG10(SOF)
  370    FORMAT('  Fitted spectral background used: Log SOF = ',F7.4)
  371    FORMAT('                                       ROH =  ',F5.3)
         WRITE(*,*)
         WRITE(*,370)XLSOF
         WRITE(*,371)ROHF
!   Having found best fit now calculate best fit AR(1) spectrum [for
!   irregularly spaced data this means out to Nyquist based on minimum
!   sample interval (instead of mean SI Nyquist used in fitting)]. 
         ROHF2=ROHF*ROHF
         DO 375 I=2,MN
           ARP2(I)=SOF*(1.0-ROHF2)/(1.0-2*ROHF*COS((FREQ(I)/XNYQ)*XPI)
     *     +ROHF2)
           XLARP2(I)=ALOG10(ARP2(I))
  375    ENDDO
       ENDIF
