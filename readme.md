
#### install dependencies

install dependencies with mamba (or conda):

`mamba env  create -f environment.yml`

#### create meshes

to generate the meshes, run: 

```
python3 scripts/extract_surfaces.py
python3 scripts/extract_vessels.py
python3 scripts/generate_mesh.py
```