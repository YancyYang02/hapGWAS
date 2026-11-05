# SACJ simulation (v0.4)

n=500, n_rep=300, LMM h2=0.3, α=0.05

| scenario | joint | ACJ | **SACJ** | alleleBon | snpBon |
|---|---:|---:|---:|---:|---:|
| null_neighbor | 0.027 | 0.043 | **0.053** | 0.027 | 0.050 |
| null_dist2 | 0.030 | 0.050 | **0.030** | 0.027 | 0.043 |
| power_single | 1.000 | 0.970 | **0.997** | 0.460 | 0.657 |
| power_neighbor_h1 | 0.923 | 0.950 | **0.870** | 0.620 | 0.740 |
| power_dist2 | 0.987 | 0.870 | **0.910** | 0.793 | 0.423 |

SACJ vs ACJ: ACJ merges Hamming≤1 components; SACJ soft Hamming affinity + spectral clusters (eigengap k).
Target regime for novelty: **power_dist2** (causal alleles at Hamming distance 2).
