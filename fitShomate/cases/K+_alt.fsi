T             # True if mu values are kcal/gmol, 0 if all are in kJ/gmol
T             # Have the H298 and S298 values
0             # mform (use 0 for full polynomial)
-252.14       # H298 (in kJ/gmol - always)
101.20        # S298 (in J/gmolK - always)
T             # Convert from Delta G formation (need to subtract elem entropy)
T             # convert from C to K
-0.2116865    # FAC =   (+ entropy_K(s) - 1/2 entropy_H2(g))T298) (kJ/gmol) ( mu = Delg0 + FAC)
12            # points
25     -67.51
50     -68.12
75     -68.73
100    -69.35
125    -69.98
150    -70.61
175    -71.24
200    -71.88
225    -72.52
250    -73.15
300    -74.40
350    -75.59
