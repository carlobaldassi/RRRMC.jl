# [Graph types](@id graphtype)

```@meta
CurrentModule = RRRMC
```

All graphs which can be used with the [sampling algorithms](@ref algorithms) belong to a type hierarchy.
At the top of the hierarchy, there is `AbstractGraph`:

```@docs
AbstractGraph
```

There are currently three abstract subclasses, which determine which sampling algorithms can be used:

```@docs
SimpleGraph
```

```@docs
DiscrGraph
```

```@docs
DoubleGraph
```


