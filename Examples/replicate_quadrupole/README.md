# replicate_quadrupole.py

Scales, tiles and randomly shifts the reference quadrupole of `../quadrupole_circle.yaml`, and writes
a DISCO input. `test_noise_quadrupole*.yaml` are example outputs.

```bash
python3 replicate_quadrupole.py --nrx 4 --nry 4 --noise-amp 6.6048 --noise-axes xy --seed 42 -o quad_4x4_noisy.yaml
```

Options and geometry: <https://lafourcadep.github.io/DISCO/quadrupoles/#replicate_quadrupolepy-large-tiled-samples>
