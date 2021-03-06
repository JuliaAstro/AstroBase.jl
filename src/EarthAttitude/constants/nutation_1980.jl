#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

const NUTATION_1980 = [
    # 1-10
    NutationLuniSolar(  0,  0,  0,  0,  1, -171996.0, -174.2, 0.0, 92025.0,    8.9, 0.0),
    NutationLuniSolar(  0,  0,  0,  0,  2,    2062.0,    0.2, 0.0,  -895.0,    0.5, 0.0),
    NutationLuniSolar( -2,  0,  2,  0,  1,      46.0,    0.0, 0.0,   -24.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0, -2,  0,  0,      11.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar( -2,  0,  2,  0,  2,      -3.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  1, -1,  0, -1,  0,      -3.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0, -2,  2, -2,  1,      -2.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0, -2,  0,  1,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2, -2,  2,  -13187.0,   -1.6, 0.0,  5736.0,   -3.1, 0.0),
    NutationLuniSolar(  0,  1,  0,  0,  0,    1426.0,   -3.4, 0.0,    54.0,   -0.1, 0.0),

    # 11-20
    NutationLuniSolar(  0,  1,  2, -2,  2,    -517.0,    1.2, 0.0,   224.0,   -0.6, 0.0),
    NutationLuniSolar(  0, -1,  2, -2,  2,     217.0,   -0.5, 0.0,   -95.0,    0.3, 0.0),
    NutationLuniSolar(  0,  0,  2, -2,  1,     129.0,    0.1, 0.0,   -70.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  0, -2,  0,      48.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2, -2,  0,     -22.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  2,  0,  0,  0,      17.0,   -0.1, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  0,  0,  1,     -15.0,    0.0, 0.0,     9.0,    0.0, 0.0),
    NutationLuniSolar(  0,  2,  2, -2,  2,     -16.0,    0.1, 0.0,     7.0,    0.0, 0.0),
    NutationLuniSolar(  0, -1,  0,  0,  1,     -12.0,    0.0, 0.0,     6.0,    0.0, 0.0),
    NutationLuniSolar( -2,  0,  0,  2,  1,      -6.0,    0.0, 0.0,     3.0,    0.0, 0.0),

    # 21-30
    NutationLuniSolar(  0, -1,  2, -2,  1,      -5.0,    0.0, 0.0,     3.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  0, -2,  1,       4.0,    0.0, 0.0,    -2.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  2, -2,  1,       4.0,    0.0, 0.0,    -2.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  0, -1,  0,      -4.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  2,  1,  0, -2,  0,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0, -2,  2,  1,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1, -2,  2,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  0,  0,  2,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  0,  1,  1,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  2, -2,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),

    # 31-40
    NutationLuniSolar(  0,  0,  2,  0,  2,   -2274.0,   -0.2, 0.0,   977.0,   -0.5, 0.0),
    NutationLuniSolar(  1,  0,  0,  0,  0,     712.0,    0.1, 0.0,    -7.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2,  0,  1,    -386.0,   -0.4, 0.0,   200.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2,  0,  2,    -301.0,    0.0, 0.0,   129.0,   -0.1, 0.0),
    NutationLuniSolar(  1,  0,  0, -2,  0,    -158.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  2,  0,  2,     123.0,    0.0, 0.0,   -53.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  0,  2,  0,      63.0,    0.0, 0.0,    -2.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  0,  0,  1,      63.0,    0.1, 0.0,   -33.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  0,  0,  1,     -58.0,   -0.1, 0.0,    32.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  2,  2,  2,     -59.0,    0.0, 0.0,    26.0,    0.0, 0.0),

    # 41-50
    NutationLuniSolar(  1,  0,  2,  0,  1,     -51.0,    0.0, 0.0,    27.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2,  2,  2,     -38.0,    0.0, 0.0,    16.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  0,  0,  0,      29.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2, -2,  2,      29.0,    0.0, 0.0,   -12.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  2,  0,  2,     -31.0,    0.0, 0.0,    13.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2,  0,  0,      26.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  2,  0,  1,      21.0,    0.0, 0.0,   -10.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  0,  2,  1,      16.0,    0.0, 0.0,    -8.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  0, -2,  1,     -13.0,    0.0, 0.0,     7.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  2,  2,  1,     -10.0,    0.0, 0.0,     5.0,    0.0, 0.0),

    # 51-60
    NutationLuniSolar(  1,  1,  0, -2,  0,      -7.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  2,  0,  2,       7.0,    0.0, 0.0,    -3.0,    0.0, 0.0),
    NutationLuniSolar(  0, -1,  2,  0,  2,      -7.0,    0.0, 0.0,     3.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2,  2,  2,      -8.0,    0.0, 0.0,     3.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  0,  2,  0,       6.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  2, -2,  2,       6.0,    0.0, 0.0,    -3.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  0,  2,  1,      -6.0,    0.0, 0.0,     3.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2,  2,  1,      -7.0,    0.0, 0.0,     3.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2, -2,  1,       6.0,    0.0, 0.0,    -3.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  0, -2,  1,      -5.0,    0.0, 0.0,     3.0,    0.0, 0.0),

    # 61-70
    NutationLuniSolar(  1, -1,  0,  0,  0,       5.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  2,  0,  1,      -5.0,    0.0, 0.0,     3.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  0, -2,  0,      -4.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0, -2,  0,  0,       4.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  0,  1,  0,      -4.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  1,  0,  0,  0,      -3.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2,  0,  0,       3.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1, -1,  2,  0,  2,      -3.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar( -1, -1,  2,  2,  2,      -3.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar( -2,  0,  0,  0,  1,      -2.0,    0.0, 0.0,     1.0,    0.0, 0.0),

    # 71-80
    NutationLuniSolar(  3,  0,  2,  0,  2,      -3.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  0, -1,  2,  2,  2,      -3.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  1,  1,  2,  0,  2,       2.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  2, -2,  1,      -2.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  0,  0,  1,       2.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  0,  0,  2,      -2.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  3,  0,  0,  0,  0,       2.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2,  1,  2,       2.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  0,  0,  2,       1.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  0, -4,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),

    # 81-90
    NutationLuniSolar( -2,  0,  2,  2,  2,       1.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  2,  4,  2,      -2.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  0, -4,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  1,  2, -2,  2,       1.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2,  2,  1,      -1.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar( -2,  0,  2,  4,  2,      -1.0,    0.0, 0.0,     1.0,    0.0, 0.0),
    NutationLuniSolar( -1,  0,  4,  0,  2,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1, -1,  0, -2,  0,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  2, -2,  1,       1.0,    0.0, 0.0,    -1.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  2,  2,  2,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),

    # 91-100
    NutationLuniSolar(  1,  0,  0,  2,  1,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  4, -2,  2,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  3,  0,  2, -2,  2,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0,  2, -2,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  2,  0,  1,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar( -1, -1,  0,  2,  1,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0, -2,  0,  1,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2, -1,  2,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  0,  2,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0, -2, -2,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),

    # 101-106
    NutationLuniSolar(  0, -1,  2,  0,  1,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  1,  0, -2,  1,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  1,  0, -2,  2,  0,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  2,  0,  0,  2,  0,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  0,  2,  4,  2,      -1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
    NutationLuniSolar(  0,  1,  0,  1,  0,       1.0,    0.0, 0.0,     0.0,    0.0, 0.0),
]

