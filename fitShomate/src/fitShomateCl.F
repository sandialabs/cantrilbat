*
*
*
* $Id: fitShomateCl.F,v 1.1 2006/01/12 00:14:47 hkmoffa Exp $

      PROGRAM fitshomate
c Maximum number of least squares conditions to be satisfied
      integer    IP1 
      parameter (IP1 = 200)
      INTEGER ii
      DOUBLE PRECISION T(IP1), mu0(IP1)
      DOUBLE PRECISION m_coeff(7), H298, S298, mu298, fac
      DOUBLE PRECISION Diff, DiffM, mu, Cp, H, S, T298
      INTEGER mform, i, Imax
      LOGICAL have298
      DOUBLE PRECISION Tmin, Tmax

      mform = 0
      H298 = 0.0
      S298 = 0.0

      ii = 12
      T(1) = 25.
      T(2) = 50.
      T(3) = 75.
      T(4) = 100.
      T(5) = 125.
      T(6) = 150.
      T(7) = 175.
      T(8) = 200.
      T(9) = 225.
      T(10) = 250.
      T(11) = 300.
      T(12) = 350.
      DO i = 1, ii
         T(i) = T(i) + 273.15
      END DO
*
* Values in kcal/gmol from Supcrt92 J PhysChemRef
      mu0(1) = -31.38
      mu0(2) = -31.69
      mu0(3) = -31.95
      mu0(4) = -32.17
      mu0(5) = -32.33
      mu0(6) = -32.45
      mu0(7) = -32.51
      mu0(8) = -32.51
      mu0(9) = -32.44
      mu0(10) = -32.27
      mu0(11) = -31.51
      mu0(12) = -29.54
*
* convert from kcal to kJ
*
      DO i = 1, ii
         mu0(i) = mu0(i) * 4.184
      END DO
*
* convert from Gibbs free energies of formation to absolute
* Gibbs free energies (NIST basis)
*
      FAC = 0.5 * (38.96224 + 66.5116)
      DO i = 1, ii
         mu0(i) = mu0(i) - FAC
      END DO
      DO i = 1, ii
         print *, T(i) , mu0(i)
      END DO
      H298 = -167.08
      S298 = 56.60
      mu298 = H298 - S298 * 298.15D0 / 1000.D0
      print *, 298.15, mu298
      mu0(1) = mu298
      have298 = .TRUE.
      mform = 0

      Tmax = T(1)
      Tmin = T(1)
      DO i = 2, ii
         IF (Tmax .LT. T(i)) THEN
            Tmax = T(i)
         ENDIF
         IF (Tmin .GT. T(i)) THEN
            Tmin  = T(i)
         ENDIF
      END DO
      IF (have298) THEN
        IF (Tmin .GT. 298.0) THEN
           Tmin = 298.0
        ENDIF
      ENDIF

      call fitshomatepoly(have298, mform, H298, S298, ii, T, 
     1                    mu0, m_coeff)

*
*          Check for Satisfaction of S and H values at T = 298.15
*
      T298 = 298.15
      CALL SH_THERM (T298, m_coeff, Cp, H, S, mu)
      
      PRINT 9995, H298, H
      PRINT 9996, S298, S
*
*  Check for Agreement with G - find the maximum deviation
*
      DiffM = 0.0
      DO I = 1, ii
        CALL SH_THERM (T(i), m_coeff, Cp, H, S, mu) 
        Diff = ABS(mu - mu0(i))
        IF (Diff .GT. DiffM) THEN
          Imax = i
          DiffM = Diff
        END IF
        print 9997, T(i), mu0(i), mu
      END DO
*
      PRINT *
      CALL SH_THERM (T(Imax), m_coeff, Cp, H, S, mu) 
      print 9988, DiffM, T(Imax), mu0(Imax), mu
      PRINT *


      PRINT *
      PRINT *,'  <thermo>'
      print 4001, Tmax, Tmin
 4001 FORMAT( '     <Shomate Pref="1 bar" Tmax="',F9.2,
     1        '" Tmin="', F9.2, '">')
      PRINT *,'        <floatArray size="7">'
      PRINT 4002, m_coeff(1), m_coeff(2), m_coeff(3)
 4002 FORMAT( '          ', 1PG17.10, ', ', 1PG17.10, ', ',
     1       1PG17.10, ', ')
      PRINT 4002, m_coeff(4), m_coeff(5), m_coeff(6)
      PRINT 4003, m_coeff(7)
 4003 FORMAT( '          ', 1PG17.10)
      PRINT *,'         </floatArray>'
      PRINT *,'      </Shomate>'
      PRINT *,'  </thermo>'
*
      RETURN
********************************************************************
9990  format (//' T         mu_data          mu_calc')
 9997 format (f10.3, 1PG15.5, 1PG15.5)

9995  format (//'Check for Satisfaction of H298 condition: '/
     $        10X, 'H298_input = ', 1PG20.13, ' cal/mole'/
     $        10X, 'H298_calc  = ', 1PG20.13, ' cal/mole')
9996  format (//'Check for Satisfaction of S298 condition: '/
     $        10X, 'S298_input = ', 1PG20.13, ' cal/moleK'/
     $        10X, 'S298_calc  = ', 1PG20.13, ' cal/moleK')

 9988 FORMAT('Max Dev: ', 1PG10.4, 'kJ/gmol T = ',  1PG10.4,
     $       ': ', 1PG15.5, ', ', 1PG15.5)
*
      END
C
C ----------------------------------------------------------------------
C
