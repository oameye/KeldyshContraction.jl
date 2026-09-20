using Preferences: set_preferences!
set_preferences!("KeldyshContraction", "dispatch_doctor_mode" => "error")

using KeldyshContraction
using Test

include("physical_twopi_response_aware_wigner.jl")
include("physical_twopi_response_aware_collision.jl")
include("physical_twopi_response_aware_weak_oracle.jl")
