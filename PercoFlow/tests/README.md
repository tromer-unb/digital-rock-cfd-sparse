# tests/

```bash
pytest -q                     # everything (about 1.5 min on one core, OpenFOAM included)
pytest -q -m "not openfoam"   # finite-element part only
pytest -q -m openfoam         # OpenFOAM cross-check only
pytest -q -m regression       # reference numbers of the example inputs
```

`conftest.py` limits BLAS/OpenMP threads to 2, so that the suite does not occupy a shared machine.

## What each test protects

| file | test | what would be wrong if it failed |
|---|---|---|
| `test_core.py` | `test_config_roundtrip_and_typos` | options lost when written and read back, or a typo accepted silently |
| | `test_band_and_crop` | wrong analysed window |
| | `test_percolation_classification` | percolating, dead-end and isolated clusters mixed up |
| | `test_critical_throat_of_a_channel` | wrong critical throat width |
| | `test_contours_square_with_hole` | solid inclusions not treated as holes |
| | `test_reservoirs` | wrong reservoir width or position |
| | `test_poiseuille_exact_through_whole_pipeline` | any error in assembly, boundary conditions, flux integration or pressure probing (Taylor–Hood must be exact); energy identity broken; Mini not below the exact value |
| | `test_tb_equals_rotated_lr` | the result depends on the flow-direction convention |
| | `test_minres_true_residual_gate` | MINRES accepts a solution whose true residual is too large |
| | `test_cylinder_array_vs_hasimoto` | the finite-element solver disagrees with the analytic series |
| `test_openfoam.py` | `test_openfoam_is_installed` | OpenFOAM missing (the test fails with instructions; it is never skipped) |
| | `test_face_connectivity_rule` | corner contacts passed to the finite-volume mesh |
| | `test_channel_exact` | the OpenFOAM pixel route or the probe-line post-processing is wrong |
| | `test_image_pipeline_with_crosscheck` | the image run, its cross-check or its manifest is broken |
| | `test_mesh_pipeline_with_crosscheck` | the mesh run or its cross-check is broken |
| | `test_disabling_is_recorded` | the cross-check can be switched off without a trace |
| `test_regression.py` | `test_example_image_is_the_published_one` | the example generator no longer reproduces the stored image |
| | `test_rock_image_reference_values` | a code change modified the computed numbers of the image route |
| | `test_obstacle_mesh_reference_values` | a code change modified the computed numbers of the mesh route |
