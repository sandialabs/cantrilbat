*
*
*
* $Id: fitShomate.F,v 1.3 2008/10/17 19:28:42 hkmoffa Exp $

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
      logical convertKcal
      logical convertFromGformation, convertCtoK


*
* Read the input file in a very rigid form
      read(5,*) convertKcal
      read(5,*) have298
      mform = 0
      read(5,*) mform
      read(5,*) H298
      read(5,*) S298
      read(5,*) convertFromGformation
      read(5,*) convertCtoK
      read(5,*) FAC
      read(5,*) ii
     
      DO i = 1, II
         read(5,*) T(i), mu0(i)
      END DO
      IF (convertCtoK) THEN
         DO i = 1, ii
            T(i) = T(i) + 273.15
         END DO
      ENDIF
*
* convert from kcal to kJ
*
      IF (convertKcal) THEN
         DO i = 1, ii
            mu0(i) = mu0(i) * 4.184
         END DO
      ENDIF
*
* convert from Gibbs free energies of formation to absolute
* Gibbs free energies (NIST basis)
*
*      FAC = - 0.5 * (38.96224 + 66.5116)
      IF (convertFromGformation) THEN
         DO i = 1, ii
            mu0(i) = mu0(i) + FAC
         END DO
      ENDIF
      DO i = 1, ii
         print *, T(i) , mu0(i)
      END DO

      IF (have298) THEN
         mu298 = H298 - S298 * 298.15D0 / 1000.D0
         print *, 298.15, mu298
         IF (ABS(T(1) - 298.15) .LT. .1) THEN
            mu0(1) = mu298
         ENDIF
      ENDIF
 

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
