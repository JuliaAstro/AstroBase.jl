#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

using ..Time: Epoch, Period, unit, value
using ..Util: CubicSpline

export TimeSeries

struct TimeSeries{S,ET,PT<:Period,T} <: AbstractArray{T,1}
    epoch::Epoch{S,ET}
    time::Vector{PT}
    data::Vector{T}
    interp::CubicSpline{Float64,T}
    function TimeSeries(ep::Epoch{S,ET}, t, u) where {S,ET}
        time_array = collect(t)
        PT = eltype(time_array)
        T = eltype(u)
        interp = CubicSpline(float(value.(time_array)), u)
        new{S,ET,PT,T}(ep, time_array, u, interp)
    end
end

function (ts::TimeSeries)(p::Period)
    p < ts.time[1] && throw(ArgumentError("`p` is too small"))
    p > ts.time[end] && throw(ArgumentError("`p` is too large"))
    return ts.interp(value(p))
end

function (ts::TimeSeries)(ep::Epoch)
    ep < ts.epoch && throw(ArgumentError("`ep` is too small"))
    ep > (ts.epoch + ts.time[end]) && throw(ArgumentError("`ep` is too large"))
    return ts.interp(value(unit(ts.time[1])(ep - ts.epoch)))
end

Base.getindex(ts::TimeSeries, idx) = ts.data[idx]
Base.size(ts::TimeSeries) = size(ts.data)
Base.eltype(ts::TimeSeries) = eltype(ts.data)
