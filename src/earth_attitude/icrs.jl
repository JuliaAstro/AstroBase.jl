export cip_coords

"""
    cip_coords(rbpn)

Extract from the bias-precession-nutation matrix the X,Y coordinates
of the Celestial Intermediate Pole.

# References

- [ERFA](https://github.com/liberfa/erfa/blob/master/src/bpn2xy.c)
"""
cip_coords(rbpn) = rbpn[3,1], rbpn[3,2]

