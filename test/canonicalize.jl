using KeldyshContraction
using KeldyshContraction: canonicalize, Bulk, In, Out

@qfields c::Destroy(Classical) q::Destroy(Quantum)

# Test the example from your request
vs5 = [
    (c(Out()), q'(Bulk(1))),
    (c(Bulk(1)), q'(Bulk(2))),
    (c(Bulk(1)), q'(Bulk(3))),
    (c(Bulk(3)), q'(Bulk(1))),
    (c(Bulk(2)), q'(Bulk(3))),
    (c(Bulk(3)), q'(Bulk(2))),
    (c(Bulk(2)), q'(In())),
]

vs6 = [
    (c(Out()), q'(Bulk(1))),
    (c(Bulk(1)), q'(Bulk(3))),
    (c(Bulk(1)), q'(Bulk(2))),
    (c(Bulk(2)), q'(Bulk(1))),
    (c(Bulk(3)), q'(Bulk(2))),
    (c(Bulk(2)), q'(Bulk(3))),
    (c(Bulk(3)), q'(In())),
]

canon5 = canonicalize(vs5)
canon6 = canonicalize(vs6)
@test canon5 == canon6

all_positions = vcat([collect(KeldyshContraction.position.(c)) for c in vs5]...)
