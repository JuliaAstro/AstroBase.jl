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

