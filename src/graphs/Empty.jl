module Empty

using ExtractMacro
using ..Interface
using ..Common
using ..QT

export GraphEmpty

import ..Interface: energy, delta_energy

type GraphEmpty <: SimpleGraph{Int}
    N::Int
    GraphEmpty(N::Integer) = return new(N)
end

energy(X::GraphEmpty, C::Config) = 0
delta_energy(X::GraphEmpty, C::Config, move::Int) = 0

end # module
