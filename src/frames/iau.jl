using ..Bodies
using ..Time: TDBEpoch

using ReferenceFrameRotations: angle_to_dcm, ddcm

for name in ALL_NAMES
    endswith(name, "Barycenter") && continue

    frame = Symbol("IAU", name)
    cname = Symbol("iau_", lowercase(name))
    add_edge!(FRAMES, frame, :ICRF)

    @eval begin
        struct $frame <: RotatingFrame end
        const $cname = $frame()
        from_sym(::Val{$(Meta.quot(frame))}) = $cname
        export $frame, $cname

        function Rotation(::ICRF, ::$frame, ep::Epoch)
            tdb = TDBEpoch(ep)
            body = $(Symbol(name))()
            euler = euler_angles(body, tdb)
            euler′ = collect(euler_rates(body, tdb))

            m = angle_to_dcm(euler..., :ZXZ)
            m′ = ddcm(m, euler′)

            Rotation{icrf, $cname}(m, m′)
        end

        Rotation(::$frame, ::ICRF, ep::Epoch) = inv(Rotation(icrf, $cname, ep))
    end
end

