export AstroDynError

struct AstroDynError <: Exception
    msg::String
end

Base.showerror(io::IO, e::AstroDynError) = print(io, e.msg)
