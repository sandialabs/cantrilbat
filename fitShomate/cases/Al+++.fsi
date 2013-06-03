T             # True if mu values are kcal/gmol, 0 if all are in kJ/gmol
N             # Have the H298 and S298 values
0             # mform (use 0 for full polynomial)
-538.4        # H298 (in kJ/gmol - always)
-325.0        # S298 (in J/gmolK - always)
T             # Convert from Delta G formation (need to subtract elem entropy)
T             # convert from C to K
50.01466      # FAC =   (+ entropy_Al(s) - 3/2 entropy_H2(g))T298) (kJ/gmol) ( mu = Delg0 + FAC)
12            # points
25     -116.54
50     -114.49
75     -112.39
100    -110.23
125    -108.02
150    -105.75
175    -103.41
200    -100.99
225    -98.48
250    -95.86
300    -90.28
350    -85.14
