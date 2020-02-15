# [Built-in graphs](@id builtin)

```@meta
CurrentModule = RRRMC
```

Following is the list of the graph models which are provided with the module. After loading the `RRRMC` module,
they can be constructed like in this example:

```text
julia> X = RRRMC.GraphRRG(10, 3)
```

Note that for models which involve randomness in the constructor you may want to set the random seed with `Random.seed!`
before calling the constructor, for reproducibility purposes.

## Basic spin glass models

### Random regular graphs

```@docs
GraphRRG
```

```@docs
GraphRRGNormal
```

```@docs
GraphRRGNormalDiscretized
```

### Edwards-Anderson graphs

```@docs
GraphEA
```

```@docs
GraphEANormal
```

```@docs
GraphEANormalDiscretized
```

### p-spin

```@docs
GraphPSpin3
```

### Sherrington-Kirkpatrick graphs

```@docs
GraphSK
```

```@docs
GraphSKNormal
```

### Binary Neural Networks

```@docs
GraphPercStep
```

```@docs
GraphPercLinear
```

```@docs
GraphPercXEntr
```

```@docs
GraphCommStep
```

```@docs
GraphCommReLU
```

### SAT

```@docs
GraphSAT
```

## Quantum models with transverse fields

```@docs
GraphQuant
```

## Robust Ensemble models

```@docs
GraphRobustEnsemble
```

## Local Entropy models

```@docs
GraphLocalEntropy
```

```@docs
GraphTopologicalLocalEntropy
```

## Mixed models

```@docs
GraphMixed
```

## Trivial models used for testing and debugging

```@docs
GraphEmpty
```

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
GraphFieldsNormalDiscretized
```

```@docs
GraphIsing1D
```
