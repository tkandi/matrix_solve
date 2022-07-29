set -x

bsub -b -I -q q_cpc -n 6 -cgsp 64 -host_stack 1024 -share_size 15000 -ldm_share_mode 5 -ldm_share_size 128 -cache_size 32 ./matrix_solve
