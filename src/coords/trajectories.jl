#
# Copyright (c) 2018-2020 Helge Eichhorn and the AstroBase.jl contributors
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#

using Roots: find_zero

using ..Time: Period, Epoch, TimeScale

import ..Time: unit
import ..Interfaces: AbstractTrajectory

export Trajectory, initial, final, state, events, times,
    LogEntry, count_id, id, epoch, detector,
    Apocenter, Pericenter, Eclipse

function transpose_aoa(input::Vector{Vector{T}}) where T
    n = length(input[1])
    output = Vector{Vector{T}}(undef, n)
    matrix = hcat(input...)
    for i = 1:n
        output[i] = matrix[i, :]
    end
    return output
end

struct Trajectory{S,F,B,T,TT,U,PT,NT} <: AbstractTrajectory{S,F,B,T,TT}
    frame::F
    body::B
    series::NT
    events::Vector{Event{S,TT}}
    data::Vector{Vector{T}}
    time::Vector{Period{U,PT}}
    function Trajectory(epoch::Epoch{S,TT},
                        time::Vector{Period{U,PT}},
                        data::Vector{Vector{T}};
                        frame::F=icrf,
                        body::B=earth,
                        names::Vector{Symbol}=Symbol[]) where {S,F,B,T,TT,U,PT}
        columns = transpose_aoa(data)
        if isempty(names)
            n = length(columns)
            if n == 3
                names = [:x, :y, :z]
            elseif n == 6
                names = [:x, :y, :z, :vx, :vy, :vz]
            else
                throw(ArgumentError("Explicit `names` need to provided if `data` does not have either 3 or 6 columns."))
            end
        end
        series = (; zip(names, [TimeSeries(epoch, time, c) for c in columns])...)
        NT = typeof(series)
        events = Vector{Event{S, TType}}[]
        new{S,F,B,T,TT,U,PT,NT}(frame, body, series, events, data, time)
    end
end

unit(tra::Trajectory) = unit(tra.time[1])

function find_events!(tra::Trajectory, detectors::Vector{<:Detector})
    for detector in detectors
        for (i, t) in enumerate(tra.time[1:end-1])
            t1 = tra.time[i+1]
            s0 = sign(detect(detector, t, tra))
            s1 = sign(detect(detector, t1, tra))
            detected, typ = detect(s0, s1, detector)
            detected || continue

            t_event = find_zero(t -> detect(detector, t * unit(tra), tra), (value(t), value(t1)))
            evt = Event(typ, epoch(tra) + t_event * unit(tra), detector)
            evt in tra.events || push!(tra.events, evt)
        end
    end
    sort!(tra.events, lt=(a, b)->isless(a.epoch, b.epoch))
end

function detect(s0, s1, detector::DiscreteDetector)
    return s1 > s0, :discrete
end

function detect(s0, s1, detector::IntervalDetector)
    s1 > s0 && return true, :start
    s1 < s0 && return true, :stop
    return false, :none
end

function (tra::Trajectory)(t)
    typ = promote_type(eltype.(values(tra.series))...)
    out = Vector{typ}(undef, length(tra.series))
    for (i, ts) in enumerate(tra.series)
        out[i] = ts(t)
    end
    return out
end

epoch(tra::Trajectory) = tra.series[1].epoch

function keplerian(tra::Trajectory, t)
    rv = tra(t)[1:6]
    return keplerian(rv[1:3], rv[4:6], grav_param(centralbody(tra)))
end

function State(tra::Trajectory, t::Period)
    rv = tra(t)[1:6]
    return State(epoch(tra) + t, rv[1:3], rv[4:6],
                 scale=timescale(tra),
                 frame=refframe(tra),
                 body=centralbody(tra))
end

function State(tra::Trajectory, ep::Epoch)
    rv = tra(t)[1:6]
    return State(ep, rv[1:3], rv[4:6],
                 scale=timescale(tra),
                 frame=refframe(tra),
                 body=centralbody(tra))
end

function KeplerianState(tra::Trajectory, t::Period)
    ele = keplerian(tra, t)
    return KeplerianState(epoch(tra) + t, ele...,
                          scale=timescale(tra),
                          frame=refframe(tra),
                          body=centralbody(tra))
end

function KeplerianState(tra::Trajectory, ep::Epoch)
    ele = keplerian(tra, ep)
    return KeplerianState(ep, ele...,
                          scale=timescale(tra),
                          frame=refframe(tra),
                          body=centralbody(tra))
end

