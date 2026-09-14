"""
    keldysh_component(x)

Return the Keldysh component of a public propagator, self-energy, or kinetic self-energy result.
"""
keldysh_component(x::DressedPropagator) = x.keldysh
keldysh_component(x::SelfEnergy) = x.keldysh
keldysh_component(x::FourierDressedPropagator) = x.keldysh
keldysh_component(x::FourierSelfEnergy) = x.keldysh
keldysh_component(x::WignerDressedPropagator) = x.keldysh
keldysh_component(x::WignerSelfEnergy) = x.keldysh
keldysh_component(x::KineticSelfEnergy) = x.keldysh

"""
    retarded_component(x)

Return the retarded component of a public propagator, self-energy, or kinetic self-energy result.
"""
retarded_component(x::DressedPropagator) = x.retarded
retarded_component(x::SelfEnergy) = x.retarded
retarded_component(x::FourierDressedPropagator) = x.retarded
retarded_component(x::FourierSelfEnergy) = x.retarded
retarded_component(x::WignerDressedPropagator) = x.retarded
retarded_component(x::WignerSelfEnergy) = x.retarded
retarded_component(x::KineticSelfEnergy) = x.retarded

"""
    advanced_component(x)

Return the advanced component of a public propagator, self-energy, or kinetic self-energy result.
"""
advanced_component(x::DressedPropagator) = x.advanced
advanced_component(x::SelfEnergy) = x.advanced
advanced_component(x::FourierDressedPropagator) = x.advanced
advanced_component(x::FourierSelfEnergy) = x.advanced
advanced_component(x::WignerDressedPropagator) = x.advanced
advanced_component(x::WignerSelfEnergy) = x.advanced
advanced_component(x::KineticSelfEnergy) = x.advanced
