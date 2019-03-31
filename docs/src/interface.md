# [Graphs interface](@id interface)

```@meta
CurrentModule = RRRMC
```

This page contains all the functions which are needed when implementing a [graph type](@ref graphtype).
See the [built-in graphs](@ref builtin) for concrete examples (in particular, the SK, RRG and EA family of
graphs have the most complete implementations). See also the documentation for the [`Config`](@ref) type.

## Functions used by all graph types

```@docs
energy
```

```@docs
delta_energy
```

```@docs
update_cache!
```

```@docs
getN
```

```@docs
neighbors
```

## Functions used by `DiscrGraph` models

```@docs
allΔE
```

## Functions used by `DoubleGraph` models

```@docs
inner_graph
```

```@docs
delta_energy_residual
```

```@docs
update_cache_residual!
```

## Functions specific to quantum models

```@docs
Qenergy
```

```@docs
transverse_mag
```

## Functions specific to robust-ensemble models

```@docs
REenergies
```

## Functions specific to local-entropy models

```@docs
LEenergies
```

```@docs
TLEenergies
```

```@docs
cenergy
```

```@docs
distances
```
