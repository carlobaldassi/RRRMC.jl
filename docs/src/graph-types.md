# [Graph types](@id graphtype)

```@meta
CurrentModule = RRRMC
```

All graphs which can be used with the [sampling algorithms](@ref algorithms) belong to a type hierarchy.
At the top of the hierarchy, there is `AbstractGraph`:

```@docs
AbstractGraph
```

There are currently three abstract subclasses and an alias, which determine how the sampling algorithms can be used:

```@docs
SimpleGraph
```

```@docs
DiscrGraph
```

```@docs
SingleGraph
```

```@docs
DoubleGraph
```
