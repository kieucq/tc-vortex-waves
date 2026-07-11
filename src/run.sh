#!/bin/bash -l
#SBATCH -N 1
#SBATCH -J vortex-wave
#SBATCH -p gpu-debug --gpus 2
#SBATCH -t 00:59:00
##SBATCH -p gpu --gpus 4
##SBATCH -t 0:59:00
#SBATCH -A r00043
#SBATCH --mem=256G
conda deactivate
module load python/gpu/3.12.5
cd /N/u/ckieu/BigRed200/codex/vortex-waves/
set -x

nvidia-smi

torchrun --standalone --nproc-per-node=2 src/tc_gravity_waves.py --backend torch --device cuda --dtype float32 --ddp
