conda deactivate
module load python/gpu/3.12.5

torchrun --standalone --nproc-per-node=2 src/tc_gravity_waves.py --backend torch --device cuda --dtype float32 --ddp
