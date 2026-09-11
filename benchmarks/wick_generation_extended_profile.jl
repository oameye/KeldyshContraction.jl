using Profile

include(joinpath(@__DIR__, "wick_generation_profile.jl"))

println("# extended scaling")

profile_interaction("boson_g4", L_g, 4, 9)
let f = () -> DressedPropagator(L_g, Val(4), Val(9))
    f() # compile/warm up
    measurement = @timed f()
    println(
        "TIMING\tboson_g4\tseconds=",
        measurement.time,
        "\tbytes=",
        measurement.bytes,
        "\tfinal_diagrams=",
        sum(length, (measurement.value.keldysh, measurement.value.retarded, measurement.value.advanced)),
    )
end

profile_interaction("fermion_quartic3", L_f, 3, 7; simplify=false)
profile_interaction("fermion_derivative2", L_p, 2, 5; simplify=false)

println("# production flat profile: boson_g3")
DressedPropagator(L_g, Val(3), Val(7)) # ensure compiled
Profile.init(delay=0.0005)
Profile.clear()
@profile DressedPropagator(L_g, Val(3), Val(7))
Profile.print(format=:flat, sortedby=:count, mincount=3)
