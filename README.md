# Labyrinth Model Lab

An interactive 3D model of the human right membranous labyrinth (semicircular canals, ampullae, cupulae, common crus, utricle and saccule) for teaching vestibular anatomy. It runs in the browser with no install.

**Open it:** https://jakalnz.github.io/labyrinth-model-lab/

Rotate the model, switch structures on and off, recolour them, and compare the separate anatomical pieces with a single fused "skin" or the MRI fluid envelope that houses them.

It is the model-design companion to the [VOR Simulator](https://jakalnz.github.io/vor-simulator/) ([source](https://github.com/jakalnz/vor-simulator)), which uses the same geometry to demonstrate the vestibulo-ocular reflex.

## Data

All geometry comes from the IE-Map in-vivo inner-ear atlas:

> Ahmadi S-A, Raiser TM, Rühl RM, Flanagin VL, zu Eulenburg P (2021). IE-Map: a novel in-vivo atlas and template of the human inner ear. *Scientific Reports* 11, 3293. https://doi.org/10.1038/s41598-021-82716-0

`ear-data.json` holds every piece in one shared frame (millimetres, origin at the utricle centroid, right ear). It also holds the pre-built merged skin: a voxel union of the membranous pieces, blurred and re-meshed with marching cubes.

## Rebuilding the data

`scripts/build-data.mjs` regenerates `ear-data.json` from the IE-Map dataset's VTK meshes. It needs the dataset downloaded locally, with the `IEMAP_DATA` environment variable pointing to its `data/` folder, and `three` installed (for the marching-cubes tables):

```sh
npm install three
IEMAP_DATA=/path/to/IEMap_data_v_1_0/data node scripts/build-data.mjs ear-data.json
```

## Licence

Code: MIT (see `LICENSE`). Geometry: see the IE-Map publication above.
