!//////////////////////////////////////////////////////////////////////
!////  $Id: napack.f,v 1.5 2019/08/16 17:16:25 saroun Exp $
!////
!//////////////////////////////////////////////////////////////////////
!////
!////  Selection of matrix subroutines from NAPACK library (http://www.netlib.org/napack)
!////  adapted to double precision
!//////////////////////////////////////////////////////////////////////


C      ________________________________________________________
C     |                                                        |
C     |  COMPUTE THE DETERMINANT OF A GENERAL FACTORED MATRIX  |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --KFACT'S OUTPUT                         |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         KDET,E--DETERMINANT IS KDET*10.**E (E INTEGER) |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS,ALOG10,DLOG10                |
C     |________________________________________________________|
C
C------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION KDET(E,A,MA)
C------------------------------------------------
      IMPLICIT NONE
      integer :: MA
      REAL(KIND(1.D0)) A(MA),D,F,G
      REAL(KIND(1.D0)) C
      INTEGER E,H,I,J,K,L,M,N
      D = A(1)
      IF ( ABS(D) .EQ. 1236 ) GOTO 10
      WRITE(6,*) 'ERROR: MUST FACTOR WITH KFACT BEFORE COMPUTING DETERMINANT'
      STOP
10    E = 0
      IF ( D .LT. 0. ) GOTO 70
      N = A(2)
      IF ( N .EQ. 1 ) GOTO 80
      D = 1.
      F = 2.D0**64
      G = 1./F
      H = 64
      M = N + 1
      J = 0
      K = 4
      L = 3 - M + M*N
      N = L + M
      DO 40 I = K,L,M
           J = J + 1
           IF ( A(I) .GT. J ) D = -D
           IF ( A(J+N) .GT. J ) D = -D
           D = D*A(I+J)
20         IF ( ABS(D) .LT. F ) GOTO 30
           E = E + H
           D = D*G
           GOTO 20
30         IF ( ABS(D) .GT. G ) GOTO 40
           E = E - H
           D = D*F
           GOTO 30
40    CONTINUE
      D = D*A(L+M)
      IF ( E .NE. 0 ) GOTO 50
      KDET = D
      RETURN
50    IF ( D .EQ. 0.D0 ) GOTO 90
      C = DLOG10(ABS(D)) + E*DLOG10(2.D0)
      E = C
      C = C - E
      IF ( C .LE. 0.D0 ) GOTO 60
      C = C - 1
      E = E + 1
60    F = 10.D0**C
      IF ( D .LT. 0.D0 ) F = -F
      KDET = F
      RETURN
70    KDET = 0.D0
      RETURN
80    KDET = A(5)
      RETURN
90    E = 0
      GOTO 70
      END


C
C      ________________________________________________________
C     |                                                        |
C     |     FACTOR A GENERAL MATRIX WITH COMPLETE PIVOTING     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --ARRAY CONTAINING MATRIX                |
C     |                 (LENGTH AT LEAST 2 + N(N+2))           |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --MATRIX DIMENSION                       |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         A     --FACTORED MATRIX                        |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |    PACKAGE SUBROUTINES: PACK                           |
C     |________________________________________________________|
C
C------------------------------------------------
      SUBROUTINE KFACT(A,LA,N)
C------------------------------------------------
      IMPLICIT NONE
      REAL(KIND(1.D0)) :: A(*),R,S,T
      INTEGER B,C,D,E,F,G,H,I,J,K,L,LA,M,N,O,P,Q
      IF ( N .LT. LA ) CALL PACK(A,LA,N)
      R = 0.
      O = N + 1
      P = O + 1
      L = 5 + N*P
      I = -N - 3
C     ---------------------------------------------
C     |*** INSERT PIVOT ROW AND COMPUTE 1-NORM ***|
C     ---------------------------------------------
10    L = L - O
      IF ( L .EQ. 4 ) GOTO 30
      S = 0.
      DO K = 1,N
           J = L - K
           T = A(I+J)
           A(J) = T
           S = S + ABS(T)
      enddo
      IF ( R .LT. S ) R = S
      I = I + 1
      GOTO 10
30    A(1) = 1236
      A(2) = N
      A(3) = R
      Q = 3 + N*O
      I = 5 - P
      B = 0
      K = 1
40    I = I + P
      IF ( K .EQ. N ) GOTO 120
      E = N - K
      H = I
      DO M = I,Q,O
           L = I + E
C     --------------------
C     |*** FIND PIVOT ***|
C     --------------------
           DO J = M,L
             IF ( ABS(A(J)) .GT. ABS(A(H)) ) H = J
           enddo
      enddo
      C = (H-4)/O
      D = 4 + O*C + K
      G = H - D
      H = D - I
      L = I + E
      F = I - B
C     -----------------------------
C     |*** INTERCHANGE COLUMNS ***|
C     -----------------------------
      DO J = F,L
        T = A(J)
        M = J + H
        A(J) = A(M)
        A(M) = T
      enddo
      J = I - K
      A(J) = G + K
      H = G + I
      T = A(H)
      A(H) = A(I)
      A(I) = T
      B = K
      K = K + 1
      IF ( T .EQ. 0.D0 ) GOTO 120
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      M = I + 1
      DO J = M,L
        A(J) = A(J)/T
      enddo
      F = I + E*O
80    J = K + L
      H = J + G
      T = A(H)
      A(H) = A(J)
      A(J) = T
      L = E + J
      IF ( T .EQ. 0.D0 ) GOTO 100
      H = I - J
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      M = J + 1
      DO J = M,L
        A(J) = A(J) - T*A(J+H)
      enddo
100   IF ( L .LT. F ) GOTO 80
      A(L+B) = C + 1
      GOTO 40
110   A(1) = -1236
      RETURN
120   IF ( A(I) .EQ. 0 ) GOTO 110
      RETURN
      END


C
C      ________________________________________________________
C     |                                                        |
C     |   REARRANGE THE ELEMENTS OF A REAL ARRAY SO THAT THE   |
C     |  ELEMENTS OF A SQUARE MATRIX ARE STORED SEQUENTIALLY   |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         A     --REAL ARRAY CONTAINING SQUARE MATRIX    |
C     |                                                        |
C     |         LA    --LEADING (ROW) DIMENSION OF ARRAY A     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN A        |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         A     --MATRIX PACKED AT START OF ARRAY        |
C     |________________________________________________________|
C
C------------------------------------------------
      SUBROUTINE PACK(A,LA,N)
C------------------------------------------------
      IMPLICIT NONE
      REAL(KIND(1.D0)) :: A(LA)
      INTEGER H,I,J,K,L,LA,N,O

      H = LA - N
      IF ( H .EQ. 0 ) RETURN
      IF ( H .GT. 0 ) GOTO 10
      WRITE(6,*) 'ERROR: LA ARGUMENT IN PACK MUST BE .GE. N ARGUMENT'
      STOP
10    I = 0
      K = 1
      L = N
      O = N*N
20    IF ( L .EQ. O ) RETURN
      I = I + H
      K = K + N
      L = L + N
      DO J = K,L
        A(J) = A(I+J)
      enddo
      GOTO 20
      END

C------------------------------------------------
      INTEGER*4 FUNCTION KVERT(V,LV,N,W)
C invert matrix, from http://www.netlib.org/napack
C modified from KVERT by J.S.
C returns 1 if not invertible, otherwise 0
C------------------------------------------------
C      ________________________________________________________
C     |                                                        |
C     |     INVERT A GENERAL MATRIX WITH COMPLETE PIVOTING     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --WORK ARRAY WITH AT LEAST 2N ELEMENTS   |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      IMPLICIT NONE
      INTEGER*4 LV,N,H,I,J,K,L,M,O,P,Q
c      REAL*8 V(LV,1),W(1),S,T
      REAL V(LV,N),W(2*N),S,T
      K=0
      ! KVERTD(A1,NMAX,N,WK)

      IF ( N .EQ. 1 ) GOTO 120
      O = N + 1
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      Q = L
      S = ABS(V(L,L))
      DO H = L,N
        DO I = L,N
          T = ABS(V(I,H))
          IF ( T .GT. S ) then
            P = I
            Q = H
            S = T
          endif
        enddo
      enddo
      W(N+L) = P
      W(O-L) = Q
      DO I = 1,N
        T = V(I,L)
        V(I,L) = V(I,Q)
        V(I,Q) = T
      enddo
      S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0.D0 ) GOTO 130
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.D0
      S = 1.D0/S
      DO I = 1,N
        V(I,L) = -S*V(I,L)
      enddo
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0.D0 ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .NE. 0 ) then
        DO I = 1,K
          V(I,J) = V(I,J) + T*V(I,L)
        enddo
      endif
      V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO I = M,N
        V(I,J) = V(I,J) + T*V(I,L)
      enddo
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K+N)
      DO I = 1,N
        T = V(I,L)
        V(I,L) = V(I,K)
        V(I,K) = T
      enddo
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
C     --------------------
C     |*** PIVOT ROWS ***|
C     --------------------
      DO J = 1,N
        DO I = 2,N
          P = W(I)
          H = O - I
          T = V(P,J)
          V(P,J) = V(H,J)
          V(H,J) = T
        enddo
      enddo
      KVERT=0
      RETURN

120   IF ( V(1,1) .EQ. 0.D0 ) GOTO 130
      V(1,1) = 1.D0/V(1,1)
      KVERT=0
      RETURN

130   KVERT=1
      END

	  
C------------------------------------------------
      INTEGER*4 FUNCTION KVERTD(V,LV,N,W)
C invert matrix, from http://www.netlib.org/napack
C rewritten from KVERT to real*8 by J.S.
C returns 1 if not invertible, otherwise 0
C------------------------------------------------

C      ________________________________________________________
C     |                                                        |
C     |     INVERT A GENERAL MATRIX WITH COMPLETE PIVOTING     |
C     |                                                        |
C     |    INPUT:                                              |
C     |                                                        |
C     |         V     --ARRAY CONTAINING MATRIX                |
C     |                                                        |
C     |         LV    --LEADING (ROW) DIMENSION OF ARRAY V     |
C     |                                                        |
C     |         N     --DIMENSION OF MATRIX STORED IN ARRAY V  |
C     |                                                        |
C     |         W     --WORK ARRAY WITH AT LEAST 2N ELEMENTS   |
C     |                                                        |
C     |    OUTPUT:                                             |
C     |                                                        |
C     |         V     --INVERSE                                |
C     |                                                        |
C     |    BUILTIN FUNCTIONS: ABS                              |
C     |________________________________________________________|
C
      IMPLICIT NONE
      INTEGER*4 LV,N,H,I,J,K,L,M,O,P,Q
c      REAL*8 V(LV,1),W(1),S,T
      REAL*8 V(LV,N),W(2*N),S,T
      K=0
      ! KVERTD(A1,NMAX,N,WK)

      IF ( N .EQ. 1 ) GOTO 120
      O = N + 1
      L = 0
      M = 1
10    IF ( L .EQ. N ) GOTO 90
      K = L
      L = M
      M = M + 1
C     ---------------------------------------
C     |*** FIND PIVOT AND START ROW SWAP ***|
C     ---------------------------------------
      P = L
      Q = L
      S = ABS(V(L,L))
      DO H = L,N
        DO I = L,N
          T = ABS(V(I,H))
          IF ( T .GT. S ) then
            P = I
            Q = H
            S = T
          endif
        enddo
      enddo
      W(N+L) = P
      W(O-L) = Q
      DO I = 1,N
        T = V(I,L)
        V(I,L) = V(I,Q)
        V(I,Q) = T
      enddo
      S = V(P,L)
      V(P,L) = V(L,L)
      IF ( S .EQ. 0.D0 ) GOTO 130
C     -----------------------------
C     |*** COMPUTE MULTIPLIERS ***|
C     -----------------------------
      V(L,L) = -1.D0
      S = 1.D0/S
      DO I = 1,N
        V(I,L) = -S*V(I,L)
      enddo
      J = L
50    J = J + 1
      IF ( J .GT. N ) J = 1
      IF ( J .EQ. L ) GOTO 10
      T = V(P,J)
      V(P,J) = V(L,J)
      V(L,J) = T
      IF ( T .EQ. 0.D0 ) GOTO 50
C     ------------------------------
C     |*** ELIMINATE BY COLUMNS ***|
C     ------------------------------
      IF ( K .NE. 0 ) then
        DO I = 1,K
          V(I,J) = V(I,J) + T*V(I,L)
        enddo
      endif
      V(L,J) = S*T
      IF ( M .GT. N ) GOTO 50
      DO I = M,N
        V(I,J) = V(I,J) + T*V(I,L)
      enddo
      GOTO 50
C     -----------------------
C     |*** PIVOT COLUMNS ***|
C     -----------------------
90    L = W(K+N)
      DO I = 1,N
        T = V(I,L)
        V(I,L) = V(I,K)
        V(I,K) = T
      enddo
      K = K - 1
      IF ( K .GT. 0 ) GOTO 90
C     --------------------
C     |*** PIVOT ROWS ***|
C     --------------------
      DO J = 1,N
        DO I = 2,N
          P = W(I)
          H = O - I
          T = V(P,J)
          V(P,J) = V(H,J)
          V(H,J) = T
        enddo
      enddo
      KVERTD=0
      RETURN

120   IF ( V(1,1) .EQ. 0.D0 ) GOTO 130
      V(1,1) = 1.D0/V(1,1)
      KVERTD=0
      RETURN

130   KVERTD=1
      END

