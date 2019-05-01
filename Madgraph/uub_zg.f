      SUBROUTINE SUUB_ZG(P1,ANS)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> z g  
C  
C Crossing   1 is u u~ -> z g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB,     NCROSS         
      PARAMETER (NEXTERNAL=4, NCOMB= 24, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 UUB_ZG
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA GOODHEL/THEL*.FALSE./
      DATA NTRY/0/
      DATA (NHEL(IHEL,  1),IHEL=1,4) / -1, -1, -1, -1/
      DATA (NHEL(IHEL,  2),IHEL=1,4) / -1, -1, -1,  1/
      DATA (NHEL(IHEL,  3),IHEL=1,4) / -1, -1,  0, -1/
      DATA (NHEL(IHEL,  4),IHEL=1,4) / -1, -1,  0,  1/
      DATA (NHEL(IHEL,  5),IHEL=1,4) / -1, -1,  1, -1/
      DATA (NHEL(IHEL,  6),IHEL=1,4) / -1, -1,  1,  1/
      DATA (NHEL(IHEL,  7),IHEL=1,4) / -1,  1, -1, -1/
      DATA (NHEL(IHEL,  8),IHEL=1,4) / -1,  1, -1,  1/
      DATA (NHEL(IHEL,  9),IHEL=1,4) / -1,  1,  0, -1/
      DATA (NHEL(IHEL, 10),IHEL=1,4) / -1,  1,  0,  1/
      DATA (NHEL(IHEL, 11),IHEL=1,4) / -1,  1,  1, -1/
      DATA (NHEL(IHEL, 12),IHEL=1,4) / -1,  1,  1,  1/
      DATA (NHEL(IHEL, 13),IHEL=1,4) /  1, -1, -1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,4) /  1, -1, -1,  1/
      DATA (NHEL(IHEL, 15),IHEL=1,4) /  1, -1,  0, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,4) /  1, -1,  0,  1/
      DATA (NHEL(IHEL, 17),IHEL=1,4) /  1, -1,  1, -1/
      DATA (NHEL(IHEL, 18),IHEL=1,4) /  1, -1,  1,  1/
      DATA (NHEL(IHEL, 19),IHEL=1,4) /  1,  1, -1, -1/
      DATA (NHEL(IHEL, 20),IHEL=1,4) /  1,  1, -1,  1/
      DATA (NHEL(IHEL, 21),IHEL=1,4) /  1,  1,  0, -1/
      DATA (NHEL(IHEL, 22),IHEL=1,4) /  1,  1,  0,  1/
      DATA (NHEL(IHEL, 23),IHEL=1,4) /  1,  1,  1, -1/
      DATA (NHEL(IHEL, 24),IHEL=1,4) /  1,  1,  1,  1/
      DATA (  IC(IHEL,  1),IHEL=1,4) /  1,  2,  3,  4/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  36/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
      DO IHEL=1,NCOMB
          IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
             T=UUB_ZG(P ,NHEL(1,IHEL),JC(1))            
             ANS(IPROC)=ANS(IPROC)+T
              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                  GOODHEL(IHEL,IPROC)=.TRUE.
C             WRITE(*,*) IHEL,T
              ENDIF
          ENDIF
      ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION UUB_ZG(P,NHEL,IC)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u u~ -> z g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=   2,NEIGEN=  1,NEXTERNAL=4)   
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=   6, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      INCLUDE "coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[2,1,4]                                                   
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL OXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL VXXXXX(P(0,3   ),ZMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL VXXXXX(P(0,4   ),ZERO ,NHEL(4   ),+1*IC(4   ),W(1,4   ))        
      CALL FVOXXX(W(1,2   ),W(1,3   ),GZU ,ZERO    ,ZERO    ,W(1,5   ))    
      CALL IOVXXX(W(1,1   ),W(1,5   ),W(1,4   ),GG ,AMP(1   ))             
      CALL FVIXXX(W(1,1   ),W(1,3   ),GZU ,ZERO    ,ZERO    ,W(1,6   ))    
      CALL IOVXXX(W(1,6   ),W(1,2   ),W(1,4   ),GG ,AMP(2   ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)
      UUB_ZG = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          UUB_ZG =UUB_ZG+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
