"""
    keldysh_component(x)

Return the Keldysh component of a `DressedPropagator` or `SelfEnergy`.
"""
keldysh_component(x::DressedPropagator) = x.keldysh
keldysh_component(x::SelfEnergy) = x.keldysh

"""
    retarded_component(x)

Return the retarded component of a `DressedPropagator` or `SelfEnergy`.
"""
retarded_component(x::DressedPropagator) = x.retarded
retarded_component(x::SelfEnergy) = x.retarded

"""
    advanced_component(x)

Return the advanced component of a `DressedPropagator` or `SelfEnergy`.
"""
advanced_component(x::DressedPropagator) = x.advanced
advanced_component(x::SelfEnergy) = x.advanced
