<a name="unreleased"></a>
## [Unreleased]

### Bugfix
- update benches with breaking changes
- Support creating 2D transform from 3D input

### Enh
- Benchmarking the creation of a ShapeInstance

### Maint
- remove unneeded transform_ops module
- Take the transform by reference in ShapeInstance
- more stringent requirements for the Shape trait
- Generalise operator macros use for other structs

### Refactor
- Split lib.rs into smaller files


<a name="v0.1.0"></a>
## v0.1.0 - 2019-02-03
### Bug
- Import PI into context of symmetry tests

### Bugfix
- Fix overlapping shapes

### Bugfix
- Install rustfmt using cargo

### Chore
- Split shape module up into smaller components
- standardise Debug/Display across all structs
- Fix spelling errors
- Import macros explicitly instead of use_macro
- ensure no more warnings about deprecated

### Enh
- Replace Arc/Mutex shared value with raw pointers

### Maint
- remove basis from PackedState instance


[Unreleased]: https://github.com/malramsay64/packing/compare/v0.1.0...HEAD
