#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

const E0 = [
    # 1-10
    SeriesTerm(0,  0,  0,  0,  1,  0,  0,  0, 2640.96e-6, -0.39e-6),
    SeriesTerm(0,  0,  0,  0,  2,  0,  0,  0,   63.52e-6, -0.02e-6),
    SeriesTerm(0,  0,  2, -2,  3,  0,  0,  0,   11.75e-6,  0.01e-6),
    SeriesTerm(0,  0,  2, -2,  1,  0,  0,  0,   11.21e-6,  0.01e-6),
    SeriesTerm(0,  0,  2, -2,  2,  0,  0,  0,   -4.55e-6,  0.00e-6),
    SeriesTerm(0,  0,  2,  0,  3,  0,  0,  0,    2.02e-6,  0.00e-6),
    SeriesTerm(0,  0,  2,  0,  1,  0,  0,  0,    1.98e-6,  0.00e-6),
    SeriesTerm(0,  0,  0,  0,  3,  0,  0,  0,   -1.72e-6,  0.00e-6),
    SeriesTerm(0,  1,  0,  0,  1,  0,  0,  0,   -1.41e-6, -0.01e-6),
    SeriesTerm(0,  1,  0,  0, -1,  0,  0,  0,   -1.26e-6, -0.01e-6),

    # 11-20
    SeriesTerm(1,  0,  0,  0, -1,  0,  0,  0,   -0.63e-6,  0.00e-6),
    SeriesTerm(1,  0,  0,  0,  1,  0,  0,  0,   -0.63e-6,  0.00e-6),
    SeriesTerm(0,  1,  2, -2,  3,  0,  0,  0,    0.46e-6,  0.00e-6),
    SeriesTerm(0,  1,  2, -2,  1,  0,  0,  0,    0.45e-6,  0.00e-6),
    SeriesTerm(0,  0,  4, -4,  4,  0,  0,  0,    0.36e-6,  0.00e-6),
    SeriesTerm(0,  0,  1, -1,  1, -8, 12,  0,   -0.24e-6, -0.12e-6),
    SeriesTerm(0,  0,  2,  0,  0,  0,  0,  0,    0.32e-6,  0.00e-6),
    SeriesTerm(0,  0,  2,  0,  2,  0,  0,  0,    0.28e-6,  0.00e-6),
    SeriesTerm(1,  0,  2,  0,  3,  0,  0,  0,    0.27e-6,  0.00e-6),
    SeriesTerm(1,  0,  2,  0,  1,  0,  0,  0,    0.26e-6,  0.00e-6),

    # 21-30
    SeriesTerm(0,  0,  2, -2,  0,  0,  0,  0,   -0.21e-6,  0.00e-6),
    SeriesTerm(0,  1, -2,  2, -3,  0,  0,  0,    0.19e-6,  0.00e-6),
    SeriesTerm(0,  1, -2,  2, -1,  0,  0,  0,    0.18e-6,  0.00e-6),
    SeriesTerm(0,  0,  0,  0,  0,  8,-13, -1,   -0.10e-6,  0.05e-6),
    SeriesTerm(0,  0,  0,  2,  0,  0,  0,  0,    0.15e-6,  0.00e-6),
    SeriesTerm(2,  0, -2,  0, -1,  0,  0,  0,   -0.14e-6,  0.00e-6),
    SeriesTerm(1,  0,  0, -2,  1,  0,  0,  0,    0.14e-6,  0.00e-6),
    SeriesTerm(0,  1,  2, -2,  2,  0,  0,  0,   -0.14e-6,  0.00e-6),
    SeriesTerm(1,  0,  0, -2, -1,  0,  0,  0,    0.14e-6,  0.00e-6),
    SeriesTerm(0,  0,  4, -2,  4,  0,  0,  0,    0.13e-6,  0.00e-6),

    # 31-33
    SeriesTerm(0,  0,  2, -2,  4,  0,  0,  0,   -0.11e-6,  0.00e-6),
    SeriesTerm(1,  0, -2,  0, -3,  0,  0,  0,    0.11e-6,  0.00e-6),
    SeriesTerm(1,  0, -2,  0, -1,  0,  0,  0,    0.11e-6,  0.00e-6),
]

const E1 = [
    SeriesTerm(0,  0,  0,  0,  1,  0,  0,  0,    -0.87e-6,  0.00e-6),
]

