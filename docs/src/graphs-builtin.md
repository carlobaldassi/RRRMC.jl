# [Built-in graphs](@id builtin)

```@meta
CurrentModule = RRRMC
```

Following is the list of the graph models which are provided with the module. After loading the `RRRMC` module,
they can be constructed like in this example:

```text
julia> X = RRRMC.GraphRRG(10, 3)
```

Note that for models which involve randomness in the constructor you may want to set the random seed with `srand`
before calling the constructor, for reproducibility purposes.

## Spin glass models

### Random regular graphs

```@docs
GraphRRG
```

```@docs
GraphRRGCont
```

```@docs
GraphRRGContSimple
```

### Edwards-Anderson graphs


```@docs
GraphEA
```

```@docs
GraphEACont
```

```@docs
GraphEAContSimple
```

### p-spin

```@docs
GraphPSpin3
```
### Quantum models with transverse fields

```@docs
GraphQIsingT
```

## Trivial models used for testing and debugging

```@docs
GraphTwoSpin
```

```@docs
GraphThreeSpin
```

```@docs
GraphFields
```

```@docs
GraphFieldsCont
```

```@docs
GraphIsing1D
```

```@docs
GraphQ0T
```


