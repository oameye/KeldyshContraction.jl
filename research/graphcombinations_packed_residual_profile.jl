using BenchmarkTools
import GraphCombinations as GC

include(ENV["KC_CORPUS_DEFS"])

function _single_atom(expression)
    sector, polynomial = only(KC.occupation_reduced_terms(expression))
    monomial, _ = only(polynomial)
    kinematic_monomial, _ = only(KC.kinematic_factor(sector))
    return KC._kinematic_atom_sector(sector, kinematic_monomial), monomial
end

function _profile_case(label, expression)
    graph_capacity, loop_capacity = KC._loop_gc_capacities(expression)
    atom_sector, monomial = _single_atom(expression)
    workspace = KC._LoopGCQuotientWorkspace(graph_capacity, loop_capacity)

    KC._graphcombinations_projective_canonical_loop_transform(atom_sector, monomial, workspace)
    storage = workspace.canonicalization
    KC._canonicalize_loop_builder!(storage, workspace.builder)

    capacity = median(
        @benchmark KC._loop_gc_capacities($expression) samples = 101 evals = 1
    )
    workspace_ctor = median(
        @benchmark KC._LoopGCQuotientWorkspace($graph_capacity, $loop_capacity) samples = 101 evals = 1
    )
    packed_search = median(
        @benchmark GC.canonicalize_directed_packed!(
            $(storage.buffer), $(storage.workspace), $(storage.graph), $(storage.labels)
        ) samples = 101 evals = 1
    )
    builder_canonicalization = median(
        @benchmark KC._canonicalize_loop_builder!($storage, $(workspace.builder)) samples = 101 evals = 1
    )
    preallocated_transform = median(
        @benchmark KC._graphcombinations_projective_canonical_loop_transform(
            $atom_sector, $monomial, $workspace
        ) samples = 101 evals = 1
    )
    full_gc = median(
        @benchmark KC.quotient_loop_momenta($expression) samples = 101 evals = 1
    )
    full_nauty = median(
        @benchmark nauty_quotient_loop_momenta($expression) samples = 101 evals = 1
    )

    us(x) = round(x.time / 1.0e3; digits=2)
    println(
        label,
        " | capacity=",
        graph_capacity,
        "v/",
        loop_capacity,
        "l | scan=",
        us(capacity),
        "us/",
        capacity.memory,
        "B | workspace=",
        us(workspace_ctor),
        "us/",
        workspace_ctor.memory,
        "B | packed-search=",
        us(packed_search),
        "us/",
        packed_search.memory,
        "B | builder+canon=",
        us(builder_canonicalization),
        "us/",
        builder_canonicalization.memory,
        "B | prealloc-transform=",
        us(preallocated_transform),
        "us/",
        preallocated_transform.memory,
        "B | full-GC=",
        us(full_gc),
        "us/",
        full_gc.memory,
        "B | Nauty=",
        us(full_nauty),
        "us/",
        full_nauty.memory,
        "B",
    )
end

profile_cases = (
    "boson/shell/pair/dense" => corpus_expression(Boson, gc_corpus_ϕ, 2, :shell, :pair, :dense),
    "fermion/shell/none/dense" =>
        corpus_expression(Fermion, gc_corpus_ψ, 2, :shell, :none, :dense),
    "boson/none/none/dense" => corpus_expression(Boson, gc_corpus_ϕ, 2, :none, :none, :dense),
    "boson/none/none/sparse" =>
        corpus_expression(Boson, gc_corpus_ϕ, 2, :none, :none, :sparse),
    "boson/none/pair/coupled" =>
        corpus_expression(Boson, gc_corpus_ϕ, 2, :none, :pair, :coupled),
    "boson/shell/quartic/sparse" =>
        corpus_expression(Boson, gc_corpus_ϕ, 2, :shell, :quartic, :sparse),
)

println("packed GC residual-cost profile")
for (label, expression) in profile_cases
    _profile_case(label, expression)
end
