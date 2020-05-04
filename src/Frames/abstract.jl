#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

export AbstractFrame, InertialFrame, RotatingFrame, isinertial, isrotating

"""
    AbstractFrame

Abstract supertype for all reference frames.
"""
abstract type AbstractFrame end

Base.show(io::IO, frame::AbstractFrame) = print(io, string(nameof(typeof(frame))))

"""
    InertialFrame

Abstract supertype for (pseudo-)inertial reference frames.
"""
abstract type InertialFrame <: AbstractFrame end

"""
    RotatingFrame

Abstract supertype for rotating reference frames.
"""
abstract type RotatingFrame <: AbstractFrame end

"""
    isinertial(frame)

Return `true` if `frame` is a (pseudo-)inertial reference frame.
"""
isinertial(::InertialFrame) = true
isinertial(::AbstractFrame) = false

"""
    isrotating(frame)

Return `true` if `frame` is a rotating reference frame.
"""
isrotating(::RotatingFrame) = true
isrotating(::AbstractFrame) = false

