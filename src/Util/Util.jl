#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

module Util

using LinearAlgebra: Tridiagonal, norm, ⋅, ×, normalize, dot
using StaticArrays: SVector

export days_to_hms, deg_to_dms, deg_to_rad, deg_to_sec, dms_to_deg, dms_to_rad,
    hms_to_days, normalize2pi, rad_to_deg, rad_to_dms, rad_to_sec, sec_to_deg,
    sec_to_rad
export angle, azimuth, elevation, vector_azel, plane_section, point_on_limb,
    spherical_to_cartesian
export CubicSpline

include("angles.jl")
include("interpolation.jl")
include("math.jl")

end

