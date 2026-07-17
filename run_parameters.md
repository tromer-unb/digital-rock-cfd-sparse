###run_horizontal
python digital_rock_cfd_sparse.py \
    --structures structures \
    --output results_horizontal \
    --direction lr \
    --boundary-band-fraction 0.08 \
    --closing-radius 0 \
    --h-min 4 \
    --h-max 30 \
    --distance-min 4 \
    --distance-max 60 \
    --viscosity 1.0 \
    --density 1.0 \
    --inlet-velocity 0.001 \
    --solver-tolerance 1e-8
###run_vertical

python digital_rock_cfd_sparse.py \
    --structures structures \
    --output results_vertical \
    --direction tb \
    --boundary-band-fraction 0.08 \
    --closing-radius 0 \
    --h-min 4 \
    --h-max 30 \
    --distance-min 4 \
    --distance-max 60 \
    --viscosity 1.0 \
    --density 1.0 \
    --inlet-velocity 0.001 \
    --solver-tolerance 1e-8

###run_both_directions
python digital_rock_cfd_sparse.py \
    --structures structures \
    --output results_both \
    --direction both \
    --boundary-band-fraction 0.08 \
    --closing-radius 0 \
    --h-min 4 \
    --h-max 30 \
    --distance-min 4 \
    --distance-max 60 \
    --viscosity 1.0 \
    --density 1.0 \
    --inlet-velocity 0.001 \
    --solver-tolerance 1e-8
