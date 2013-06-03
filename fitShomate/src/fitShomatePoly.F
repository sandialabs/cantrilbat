*
*
* $Id: fitShomatePoly.F,v 1.1 2006/01/12 00:14:47 hkmoffa Exp $
*
* ----------------------------------------------------------------
*
      subroutine fitshomatepoly(have298,
     &          mform, H298, S298, ii, T, mu0, m_coeff)
*
      INTEGER  ii,  mform
      double precision H298, S298, T(ii), mu0(ii), m_coeff(7)
*
* ----------------------------------------------------------------
*
*         Fits a Shomate polynomial to a set of Gibbs Free Energies
*    and the value of the S_298 and H_298.
*
* Parameters
*       have298 = boolean indicating that we have 298 enthalpies
*                 and entropies. If we do, we will include their
*                 values in the fit. If not, we will just use the
*                 Gibbs free energy values. 
*       S298 = value of the 25C entropy (J/gmolK)
*       H298 = value of the 25C enthalpy (kJ/gmol)
*
*       ii = number of Gibbs Free Energy Points
*       mu(ii) = value of the Gibbs free energies (kJ/gmol)
*       
*
*      Cp = A + Bt + Ct**2 + Dt**3 + E/t**2      (J/gmolK)
*
*      H - H298 = At + Bt**2/2 + Ct**3/3 + Dt**4 - E/t + F
*                                                 (kJ/gmol)
*
*      S = Alnt + Bt + Ct**2/2 + Dt**3/3 - E/(2t**2) + G
*                                                 (J/(gmol K))
*
*      where t = T/1000
*
*      Note, unlike the nasa polynomials the shomate polynomials
*      have units associated with them.
*
*      Formulation can be changed:
*
*       mform = 0 : Use all coefficients
*       mform = 1 : D = 0;
*       mform = 2   D = E = 0;
*       mform = 3   C = D = E = 0
*       mform = 4   B = C = D = E = 0
*
* Return
*       m_coeff(7) : The parameters in alphabetical order
*
* -----------------------------------------------------------------
*
*           Local Variables
*
      integer    LP, NCP, LDAT, I1, I2, I3, I4
      logical    have298
      parameter (LP = 7, NCP = 5,
     &           LDAT = 20, I1 = 1, I2 = 2, I3 = 3, I4 = 4)
c
c Maximum number of least squares conditions to be satisfied
      integer    IP1 
      parameter (IP1 = 200)
c Maximum number of total conditions to be satisfied
      integer    IDIM
      parameter (IDIM = IP1 + 2)
c Maximum number of parameters to be found
      integer    MAXPAR0
      parameter (MAXPAR0 = 10)
c Maximum number of parameters plus one
      integer    JDIM
      parameter (JDIM = MAXPAR0 + 1)
c Integer work space requirements
      integer    IP2
      parameter (IP2 = 2*MAXPAR0 + 2)
c Real work space reqirements
      integer    IWS0
      parameter (IWS0 = 2*2 + 4*MAXPAR0 + IP1 + 14 )
c
C Number of parameters in the least squares fit
      INTEGER npar
*
      integer i, kdate, n, j, nn, l, 
     &        me, ma, mg,  mode
  
      integer iarray(3), ip(IP2)
      DOUBLE PRECISION cmin, cmax
      DOUBLE PRECISION tr, trp2, trp3, trp4
      DOUBLE PRECISION  x(MAXPAR0), w(IDIM, JDIM), prgopt(13), 
     &        ws(IWS0)
      DOUBLE PRECISION  rnorme, rnorml
*
      DOUBLE PRECISION D1MACH
      EXTERNAL         D1MACH

C
C           EXTERNALS
C
      integer  ilast
      external ilast
*

      if (ii .gt. IP1) then
        PRINT *, 'fitShomatePoly: FATAL ERROR, ',
     &         'not enough storage space allocated'
        return
      END IF
      IF (mform .LT. 0 .OR. mform .GT. 4) THEN
         PRINT *,'fitShomatePoly: Invalid mform: ', mform
         STOP 3
      END IF 
*
*              Change the data to single precision
*

*
*            Set the program options for DLSEI fitting program
*
      prgopt(1) = 4
      prgopt(2) = 2
      prgopt(3) = 1.0
*
*            Change the tolerance on the equality constraints
*            to the minimum values, possible
*
      prgopt(4) = 7
      prgopt(5) = 4
      prgopt(6) = D1MACH(4)
*
*            Change the tolerance on the minimization problem
*
      prgopt(7) = 10
      prgopt(8) = 5
      prgopt(9) = D1MACH(4)
      prgopt(10) = 1
*
* get the date
C*****sun
      call idate(iarray)
      iarray(3) = mod(iarray(3),100)
      kdate = ((iarray(2)*10000) + iarray(1)*100 + iarray(3) )
C*****END sun
*
* Figure out the number of parameters
      IF (mform .EQ. 0) THEN
         npar = 7
      ELSE
         npar = 7 - mform
      ENDIF
*
* Figure out the number of requirements
*
      IF (have298) THEN
        n = ii + 2
      ELSE
        n = ii
      ENDIF
*
* Check to see that we are at least over or equally 
* constrained
*
      IF (have298) THEN
        IF (II .LT. 2) THEN
           PRINT *,'fitShomatePoly: II must be equal to 2 or greater'
           STOP 3
        ENDIF
      ELSE
        IF (II .LT. 3) THEN
           PRINT *,'fitShomatePoly: II must be equal to 3 or greater'
           STOP 3
        ENDIF
      ENDIF
      IF (npar .GT. n) THEN
        mform = 7 - n
        npar = 7 - mform
      ENDIF
*     
* First check for degenerate conditions
*
      CALL minmax(mu0, ii, cmin, cmax)
      IF (cmin .eq. cmax) then
        DO i = 1, 7
          m_coeff(i) = 0.0
        END DO
        m_coeff(6) = cmin
      ELSE
*
* Then do the general problem
*
        DO 95 i = 1, idim
          DO 95 j = 1, jdim
            w(i,j) = 0.0d0
   95   continue

        nn = npar + 1
        L = 0
*
* Do the first H298 equation
        IF (have298) THEN
           L = L + 1
           w(L,nn) = H298
           tr = 298.15 / 1000.
           trp2 = tr * tr
           trp3 = trp2 * tr
           trp4 = trp3 * tr
           W(L,1) = tr
           IF (mform .LT. 4) THEN
              W(L,2) = trp2 / 2.0
           ENDIF
           IF (mform .LT. 3) THEN
              W(L,3) = trp3/3.0
           ENDIF
           IF (mform .EQ. 1) THEN
              W(L,4) = - 1.0 / tr
           ELSE IF (mform .EQ. 0) THEN
              W(L,4) = trp4/ 4.0
              W(L,5) = - 1.0 / tr
           ENDIF
           W(L,npar-1) = 1.0
           W(L,npar) = 0.0
        ENDIF
*
* Do the S298 Equation
        IF (have298) THEN
           L = L + 1
           w(L,nn) = S298
           W(L,1) = log(tr)
           IF (mform .LT. 4) THEN
              W(L,2) = tr
           ENDIF
           IF (mform .LT. 3) THEN
              W(L,3) = trp2/2.0
           ENDIF
           IF (mform .EQ. 1) THEN
              W(L,4) = - 1.0 / (trp2 * 2.0)
           ELSE IF (mform .EQ. 0) THEN
              W(L,4) = trp3/ 3.0
              W(L,5) = - 1.0 / (trp2 * 2.0)
           ENDIF
           W(L,npar-1) = 0.0
           W(L,npar)   = 1.0
        END IF
*
* Do the mu vs T equations
*
        DO i = 1, ii
           L = L + 1
           w(L,nn) = mu0(i)
           tr = T(i) / 1000.
           trp2 = tr * tr
           trp3 = trp2 * tr
           trp4 = trp3 * tr

           w(L,nn) = mu0(i)
           W(L,1) = tr - tr *log(tr)
           IF (mform .LT. 4) THEN
              W(L,2) = - trp2/2.0
           ENDIF
           IF (mform .LT. 3) THEN
              W(L,3) = -trp3/6.0
           ENDIF
           IF (mform .EQ. 1) THEN
              W(L,4) = - 1.0 / (tr * 2.0)
           ELSE IF (mform .EQ. 0) THEN
              W(L,4) = -trp4/12.0
              W(L,5) = - 1.0 / (tr * 2.0)
           ENDIF
           W(L,npar-1) = 1.0
           W(L,npar)   = -tr
        END DO
c     Number of equality constraints
        me = 0
c     Number of least squares equations
        ma = L
c     Number of inequality constraints
        mg = 0
c      
        ip(1) = IWS0
        ip(2) = IP2
        CALL DLSEI(w, idim, me, ma, mg, npar, prgopt, x, rnorme, rnorml,
     &            mode, ws, ip)
        IF (mode .ne. 0) then
          PRINT 7020, mode
          stop 3
        END IF

        IF (mform .EQ. 0) THEN
           DO i = 1, npar
              m_coeff(i) = x(i)
           END DO
        ELSE IF (mform .eq. 1) THEN
           DO i = 1, 3
              m_coeff(i) = x(i)
           END DO
           m_coeff(4) = 0.0
           m_coeff(5) = x(4)
           m_coeff(6) = x(5)
           m_coeff(7) = x(6)
        ELSE IF (mform .eq. 2) THEN
           m_coeff(1) = x(1)
           m_coeff(2) = x(2)
           m_coeff(3) = x(3)
           m_coeff(4) = 0.0
           m_coeff(5) = 0.0
           m_coeff(6) = x(4)
           m_coeff(7) = x(5)
        ELSE IF (mform .eq. 3) THEN
           m_coeff(1) = x(1)
           m_coeff(2) = x(2)
           m_coeff(3) = 0.0
           m_coeff(4) = 0.0
           m_coeff(5) = 0.0
           m_coeff(6) = x(3)
           m_coeff(7) = x(4)
        ELSE IF (mform .eq. 4) THEN
           m_coeff(1) = x(1)
           m_coeff(2) = 0.0
           m_coeff(3) = 0.0
           m_coeff(4) = 0.0
           m_coeff(5) = 0.0
           m_coeff(6) = x(2)
           m_coeff(7) = x(3)
        ENDIF

      end if
c

C----------------------------------------------------------------------C
C          Print Out a Sensitive check to the fitting procedure
C----------------------------------------------------------------------C
*
      call prnchk (x, H298, S298, ii, T, mu0)
*
C--------------------------------------------------------------------C
C            CLEANUP SECTION
C--------------------------------------------------------------------C
*
      return
*
C--------------------------------------------------------------------C
C               Format STATEMENTS
C--------------------------------------------------------------------C
*
 7020 format(10x,21hERROR IN DLSEI, MODE= ,i4)
*
*     Formats for writing NASA database files
*
*
      END
*
***********************************************************************
*
      SUBROUTINE SH_THERM(Tkelvin, x, Cp, H, S, mu)
*
      DOUBLE PRECISION Tkelvin, x(7), Cp, H, S, mu
*
* ----------------------------------------------------------------------
*         This routine calculates the thermo functions, 
*         Cp, H, and S
*         given a temperature, T.  Polynomial form used is the
*         shomate polynomials. It returns dimensional quantities.
*
*        Cp, S = J/(gmol K)
*        H,  G = kJ/gmol
* ----------------------------------------------------------------------
* Local Variables
*
      DOUBLE PRECISION   T, Tp2, Tp3, Tp4
      DOUBLE PRECISION   A, B, C, D, E, F, G
*
      T = Tkelvin / 1000.0
      Tp2 = T**2
      Tp3 = Tp2 * T
      Tp4 = Tp3 * T
*
      A = x(1)
      B = X(2)
      C = x(3)
      D = x(4)
      E = x(5)
      F = x(6)
      G = x(7)    
*
      Cp = A + B*T + C*Tp2 + D*Tp3 + E/Tp2 
*
      H  = A*T + B*Tp2/2.0 + C*Tp3/3.0 + D*Tp4/4.0 - E/T + F 
*
      S  = A*log(T) + B*T + C*Tp2/2.0 + D*Tp3/3.0 - E/(tp2*2.0) + G
#ifdef DO_NOT_TAKE_OUT
      IF (ABS(Tkelvin - 298.15) .LT. 0.01) THEN
        print *, 'alog(t) = ', A*log(T)
        print *, 'bt      = ', B * T
        print *, 'ctp2/2  = ', C*Tp2/2.0
        print *, 'dtp3/3  = ', D*Tp3/3.0
        print *, '-E/tp2*2= ', - E / (tp2*2.0)
        print *, 'g       = ', G
        print *, '-----------------------------'
        print *, 'S       = ', S
      END IF
#endif
*
      mu = H - T*S
*
      RETURN
      END
*
*----------------------------------------------------------------------C
*
      subroutine minmax(v, nn, vmin, vmax)
      integer              nn
      DOUBLE PRECISION  v(nn), vmin, vmax
*
*        Finds the minimum and maximum of a vector
*
      INTEGER    N
      IF (NN .GT. 0) THEN
        vmin = v(1)
        vmax = v(1)
        DO n = 2, nn
          vmax = max(v(n),vmax)
          vmin = min(v(n),vmin)
        END DO
      ELSE
        PRINT *,'MINMAX ERROR: nn lt 1, ', NN
        STOP 2
      END IF
      return
      end

*
*----------------------------------------------------------------------*
*
      integer function ifirst(istr)
      character              istr*(*)
*
*       Finds the first non-blank character in a character string
*       Returns the position of that character from the start of
*       the string.  If no non-blank characters are found, then
*       this function returns 0.
*
      integer   i
      INTRINSIC    CHAR, LEN
      do i = 1, len(istr)
        if (istr(i:i) .ne. ' '    .and.
     &      istr(i:i) .ne. '	'       ) then 
          if (istr(i:i) .eq. CHAR(0)) then
            ifirst = 0
          else
            ifirst = i
          endif
          return
        end if
      end do
      ifirst = 0
      return
      end
*
*----------------------------------------------------------------------C
*
      integer function ilast(istr)
      character              istr*(*)
*
*       Finds the last non-blank character in a character string
*       Returns the position of that character from the start of
*       the string.  If no non-blank characters are found, then
*       this function returns 0.
*
      integer   i
*
      INTRINSIC    CHAR, LEN
*
      do i = len(istr), 1, -1
        if (istr(i:i) .ne. ' '    .and.
     &      istr(i:i) .ne. '	' .and.
     &      istr(i:i) .ne. CHAR(0)    ) then
          ilast = i
          return
        end if
      end do
      ilast = 0
      return
      end
*
*----------------------------------------------------------------------*
*
      integer function STRIP(istr)
      character              istr*(*)
*
*       Strips leading and trailing blanks from 
*       a character string.  Returns the number of non-zero 
*       characters.
*       this function returns 0.
*
      integer   i, j, ifi, ila, num
*
      integer   IFIRST, ILAST
      external  IFIRST, ILAST
*
      INTRINSIC   CHAR, LEN
*
      ifi = IFIRST(istr)
      IF (IFI .EQ. 0) THEN
        strip = 0
      ELSE 
        ila = ILAST(istr)
        num = ila - ifi + 1
        IF (IFI .GT. 1) THEN
          j = ifi
          do i = 1, num
            istr(i:i) = istr(j:j)
            j = j + 1
          end do
        ENDIF
        DO i = num + 1, LEN(istr)
          istr(i:i) = CHAR(0)
        END DO
        strip = num
      ENDIF
      return
      end
*
*----------------------------------------------------------------------*
*
      subroutine prnchk (x,  h298, s298, ii, T, mu0)
*
      INTEGER           ii
      DOUBLE PRECISION  x(10), H298, S298,
     &                  T(ii), mu0(ii)
*
* Local Variables
* 
      INTEGER         I, IMAX, LT
      DOUBLE PRECISION Cp, H, S, Diff, DiffM, mu
      CHARACTER       Cnum*4
*
      DOUBLE PRECISION Rgas
      PARAMETER       (Rgas = 1.987206D0)

      DOUBLE PRECISION T298
      PARAMETER       (T298 = 298.15D0)
*
      INTEGER  STRIP
      EXTERNAL STRIP
*
*         Check the Cp, H, and S connectivity condition between the 
*         two temperature regions
*

*
*          Check for Satisfaction of S and H values at T = 298.15
*
      CALL SH_THERM (T298, x, Cp, H, S, mu)
      
      PRINT 9995, H298, H
      PRINT 9996, S298, S
*
*  Check for Agreement with G - find the maximum deviation
*
      DiffM = 0.0
      DO I = 1, ii
        CALL SH_THERM (T(i), x, Cp, H, S, mu) 
        Diff = ABS(mu - mu0(i))
        IF (Diff .GT. DiffM) THEN
          Imax = i
          DiffM = Diff
        END IF
        print 9997, T(i), mu0(i), mu
      END DO
*
 

      WRITE(Cnum,'(I4)') Imax
      LT = STRIP(Cnum)





      RETURN
********************************************************************
 9997 format (f10.3, 1PG15.5, 1PG15.5)

9995  format (//'Check for Satisfaction of H298 condition: '/
     $        10X, 'H298_input = ', 1PG20.13, ' cal/mole'/
     $        10X, 'H298_calc  = ', 1PG20.13, ' cal/mole')
9996  format (//'Check for Satisfaction of S298 condition: '/
     $        10X, 'S298_input = ', 1PG20.13, ' cal/moleK'/
     $        10X, 'S298_calc  = ', 1PG20.13, ' cal/moleK')
*
      end
C---------------------------------------------------------------------C
