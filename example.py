import pychic

Enu = 2.5 # GeV
L = 1300. # km

# Using CHIC
ch = pychic.CHIC("antineutrino")
P = ch.compute_oscillations(Enu, L)

# Using CHICDIFF
dch = pychic.CHICDIFF("antineutrino")
P = dch.compute_oscillations(Enu, L)
dP = dch.compute_oscillations_derivatives("dm221", Enu, L)

