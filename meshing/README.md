# Tools for converting meshes outputted from SimVascular to other formats
An example of using `vtk` and `meshio` to convert a mesh with boundary
markers from SimVascular to FEniCS. You get test data from the demoproject
in the [SimVascular documentation](https://simtk.org/frs/download_confirm.php/file/5113/DemoProject.zip?group_id=930, "Demo data"). You can the execute the conversion by running:

```
python convert_vtu_and_vtp.py -ivtp [path/to/vtp/file.vtp] -ivtu [path/to/vtu/file.vtu] -o [output.xdmf]
```

replacing the `[]` with your local paths.


## Dependencies
The convertion script depend on `meshio` and `vtk`, please confer with their install pages for installation instructions, but for completeness the following line will install both libraries:

```
pip install meshio vtk
```

## Help
To see all arguments which can be passed to the scripts, please run:

```
python convert_vtu_and_vtp.py -h
python example_fenics.py -h
```
