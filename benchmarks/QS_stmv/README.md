# STMV benchmark

This benchmark test the performance of CP2K to run calculations of the electronic structure of relatively complex systems containing a million atoms. The input is based on [earlier work](https://pubs.acs.org/doi/full/10.1021/acs.jctc.6b00398) where the electronic structure of the STMV virus was simulated based on DFT and subsystem DFT. Here, instead, the xTB tight-binding method is employed. The input is realistic in its input settings and might be useful to set up similar systems.

## Properties of the benchmark

The benchmark exercises in particular the sparse matrix handling and linear scaling algorithms in CP2K. It performs 1 step of geometry optimization, so requires SCF, energy and force calculations. Some properties are computed as well. Given the xTB method, relatively small block sizes dominate in the sparse matrix multiplication.

## Typical timings and setup

A typical parallel run will require on the order of 256 nodes, to see completion of the benchmark in a reasonable time, and memory consumption. An invocation with slurm on a system with 32 threads per node (dual socket, Intel(R) Xeon(R) CPU E5-2695 v4 @ 2.10GHz, Piz Daint multi-core) could look like:
```
export OMP_NUM_THREADS=8
srun --cpu-bind=none --nodes=256 --ntasks=1024 --ntasks-per-node=4 --cpus-per-task=8 ./cp2k.psmp -i stmv_xtb.inp -o stmv_xtb.out
```
Which would need roughly 7Gb per rank (28Gb per node), and would run in in less than 4h. The timing report for  this run (based on CP2K 7.0, git:bf104a630):

```
SUBROUTINE                       CALLS  ASD         SELF TIME        TOTAL TIME
                                MAXIMUM       AVERAGE  MAXIMUM  AVERAGE  MAXIMUM
 CP2K                                 1  1.0    6.012    6.487 8256.073 8256.198
 cp_geo_opt                           1  2.0    0.041    0.092 8236.532 8236.652
 geoopt_lbfgs                         1  3.0    0.010    0.074 8236.490 8236.639
 cp_opt_gopt_step                     1  4.0    0.087    0.285 8236.314 8236.439
 cp_eval_at                           2  5.0    0.097    0.148 8223.558 8223.698
 qs_forces                            2  6.0    0.831    1.322 8223.245 8223.380
 qs_energies                          2  7.0    0.034    0.117 8095.812 8099.994
 ls_scf                               2  8.0    0.029    0.093 7976.523 7980.780
 ls_scf_main                          2  9.0    0.016    0.055 7240.555 7240.672
 dm_ls_curvy_optimization            47 10.0    0.003    0.076 6583.940 6599.528
 dbcsr_multiply_generic             737 13.0    7.232    7.446 6427.778 6579.351
 multiply_cannon                    737 14.0   41.860   44.498 5883.682 6121.346
 optimization_step                   47 11.0    0.001    0.015 4866.856 4889.982
 multiply_cannon_multrec          23584 15.0 4119.121 4434.655 4133.269 4447.990
 compute_direction_newton            16 12.0    0.812    0.870 2530.534 2536.757
 update_p_exp                        47 12.0    0.009    0.015 2336.108 2359.676
 mp_waitall_1                    200904 16.1 1041.558 1716.869 1041.558 1716.869
 transform_matrix_orth               65 11.0    0.002    0.003 1664.426 1698.501
 commutator_symm                    170 13.0    0.003    0.004 1413.495 1456.582
 purify_mcweeny_orth                 50 12.9    0.003    0.065 1312.274 1334.495
 multiply_cannon_metrocomm3       23584 15.0    0.125    0.167  294.026 1106.474
 multiply_cannon_metrocomm1       23584 15.0    0.172    0.196  596.117 1034.438
 dbcsr_new_transposed               417 13.2   32.609   35.757  703.789  824.199
 dbcsr_redistribute                 417 14.2  379.063  436.415  659.380  776.357
 multiply_cannon_metrocomm4       22847 15.0    0.156    0.232  200.626  526.612
 mp_irecv_dv                      58794 16.3  191.731  503.107  191.731  503.107
 calculate_norms                  47168 15.0  437.842  473.288  437.842  473.288
 ls_scf_init_scf                      2  9.0    0.014    0.088  461.141  461.330
 ls_scf_init_matrix_S                 2 10.0    0.001    0.030  413.455  414.477
 ls_scf_dm_to_ks                     49 10.0    0.001    0.023  381.400  401.911
 matrix_sqrt_Newton_Schulz            2 11.0    0.023    0.426  389.123  389.457
 make_m2s                          1474 14.0    7.055    7.820  283.134  324.396
 mp_alltoall_i22                    569 14.7  236.860  323.813  236.860  323.813
 make_images                       1474 15.0   36.912   39.536  268.811  310.928
 ls_scf_post                          2  9.0    0.009    0.087  274.798  279.184
 mp_sum_l                          2322 13.8  207.665  270.645  207.665  270.645
 density_matrix_trs4                  1 10.0    0.015    0.042  248.757  249.218
 qs_ks_update_qs_env                 50 11.0    0.000    0.000  208.761  218.124
 make_images_data                  1474 16.0    0.081    0.122  152.463  205.226
 hybrid_alltoall_any               1525 16.9    0.731   22.546  142.014  199.699
 matrix_ls_to_qs                     51 11.0    0.001    0.001  183.276  198.338
 mp_sum_d                          2810 12.6  148.925  197.472  148.925  197.472
 rebuild_ks_matrix                   52 11.9    0.042    0.061  194.777  195.286
 build_xtb_ks_matrix                 52 12.9    5.507    6.598  194.735  195.237
 dbcsr_complete_redistribute        101 12.5   89.804   93.327  171.256  191.421
 post_scf_homo_lumo                   2 10.0    0.000    0.001  186.075  186.911
 matrix_decluster                    51 12.0    0.095    0.238  147.731  168.361
 dbcsr_finalize                    3365 14.5    2.293    2.770  140.479  155.280
 dbcsr_frobenius_norm               712 13.0   34.195   35.252  114.124  139.443
 mp_sum_dm                          430  9.1  131.588  135.181  131.588  135.181
 mp_allgather_i34                   737 15.0   71.689  130.243   71.689  130.243
 dbcsr_merge_all                   2761 15.5   72.120   75.581  118.675  123.148
 qs_energies_init_hamiltonians        2  8.0    0.007    0.082  119.151  119.333
 dbcsr_add_d                       1749 13.0    0.006    0.008  100.754  111.917
 dbcsr_add_anytype                 1749 14.0   25.357   28.939  100.748  111.910
 build_xtb_matrices                   4  8.0   64.662   80.037  110.824  111.490
 make_images_sizes                 1474 16.0    0.006    0.036   37.330  106.622
 mp_alltoall_i44                   1474 17.0   37.324  106.617   37.324  106.617
 dbcsr_dot_sd                       345 13.0   13.894   14.777   73.443  101.347
 arnoldi_extremal                     9 11.2    0.089    0.533   98.396  100.332
 arnoldi_normal_ev                    9 12.2    0.343    0.839   98.307  100.124
 build_xtb_coulomb                   52 13.9   18.489   19.417   99.082   99.771
 calculate_dispersion_pairpot         2  9.0   57.454   70.439   92.449   92.457
 ao_charges_kp_2                     52 13.9    2.220    2.512   89.914   91.795
 mp_alltoall_d11v                  1707 14.8   89.032   90.724   89.032   90.724
 ao_charges_2                        52 14.9    4.520    4.815   87.694   89.634
 build_subspace                      39 13.3    0.645    1.234   87.475   88.719
 mp_shift_i                        8184  9.0   46.858   66.427   46.858   66.427
 setup_rec_index_2d                1474 15.0   59.164   65.956   59.164   65.956
 dbcsr_matrix_vector_mult           884 14.0    0.068    1.905   61.132   64.708
 dbcsr_data_release               70614 16.3   31.961   59.053   32.162   59.250
 mp_sum_iv                          744 14.9   42.010   51.567   42.010   51.567
 mp_sum_dv                         5097 15.6   44.088   50.582   44.088   50.582
 ls_scf_initial_guess                 2 10.0    0.000    0.015   47.606   48.384
 calculate_w_matrix                   2 10.0    0.000    0.000   44.894   45.647
 dbcsr_copy_into_existing            51 12.0   35.076   44.774   35.077   44.775
 dbcsr_matrix_vector_mult_local     884 15.0   38.363   44.350   38.371   44.359
 dbcsr_multiply_generic_mpsum_f      60 12.7    0.000    0.000   28.968   43.956
 ls_scf_store_result                  2 10.0    0.020    0.206   42.334   42.976
 matrix_qs_to_ls                     50 10.0    0.001    0.021   35.255   37.254
 matrix_cluster                      50 11.0    0.115    0.300   35.254   37.254
```

and DBCSR statistics

```

 COUNTER                                    TOTAL       BLAS       SMM       ACC
 flops     2 x     2 x     2       36142856249008       0.0%    100.0%      0.0%
 flops     2 x     2 x     8       39213980112512       0.0%    100.0%      0.0%
 flops     8 x     2 x     2       43804007223424       0.0%    100.0%      0.0%
 flops     2 x     8 x     2       44584775765312       0.0%    100.0%      0.0%
 flops     4 x     2 x     2       90117651537184       0.0%    100.0%      0.0%
 flops     2 x     4 x     2       90294481986400       0.0%    100.0%      0.0%
 flops     2 x     2 x     4       97540850606240       0.0%    100.0%      0.0%
 flops     8 x     8 x     2       97563507289344       0.0%    100.0%      0.0%
 flops     2 x     4 x     8      107600053992192       0.0%    100.0%      0.0%
 flops     4 x     2 x     8      107835358956416       0.0%    100.0%      0.0%
 flops     8 x     4 x     2      117969022215680       0.0%    100.0%      0.0%
 flops     4 x     8 x     2      119092012799616       0.0%    100.0%      0.0%
 flops     8 x     2 x     4      130156579707520       0.0%    100.0%      0.0%
 flops     2 x     8 x     4      130528054439040       0.0%    100.0%      0.0%
 flops     8 x     2 x     8      161502281254912       0.0%    100.0%      0.0%
 flops     2 x     8 x     8      166307689803776       0.0%    100.0%      0.0%
 flops     4 x     4 x     2      229852669562368       0.0%    100.0%      0.0%
 flops     2 x     4 x     4      247224891759168       0.0%    100.0%      0.0%
 flops     4 x     2 x     4      249694810664896       0.0%    100.0%      0.0%
 flops     8 x     8 x     4      305672038397440       0.0%    100.0%      0.0%
 flops     4 x     4 x     8      308761351253760       0.0%    100.0%      0.0%
 flops     8 x     4 x     4      362908283625472       0.0%    100.0%      0.0%
 flops     4 x     8 x     4      365532741502720       0.0%    100.0%      0.0%
 flops     8 x     4 x     8      455052334568448       0.0%    100.0%      0.0%
 flops     4 x     8 x     8      465166240374272       0.0%    100.0%      0.0%
 flops     4 x     4 x     4      648336487473152       0.0%    100.0%      0.0%
 flops     8 x     8 x     8    10568517651483648       0.0%    100.0%      0.0%
 flops inhomo. stacks            2554975071623610     100.0%      0.0%      0.0%
 flops total                        18.341948E+15      13.9%     86.1%      0.0%
 flops max/rank                     18.975286E+12      15.8%     84.2%      0.0%
 matmuls inhomo. stacks             1752383451045     100.0%      0.0%      0.0%
 matmuls total                     55112132774528       3.2%     96.8%      0.0%
 number of processed stacks           61576905441       3.4%     96.6%      0.0%
 average stack size                                   830.2     897.3       0.0
 marketing flops                    36.676667E+21
 -------------------------------------------------------------------------------
 # multiplications                            737
 max memory usage/rank               7.564362E+09
 # max total images/rank                        1
 # max 3D layers                                1
 # MPI messages exchanged                46790656
 MPI messages size (bytes):
  total size                         4.348855E+15
  min size                           0.000000E+00
  max size                         202.967424E+06
  average size                      92.942816E+06
 MPI breakdown and total messages size (bytes):
             size <=      128              153760                        0
       128 < size <=     8192               83948                332481200
      8192 < size <=    32768               54994               1018984136
     32768 < size <=   131072              468286              34129603760
    131072 < size <=  4194304             1278316             958318296392
   4194304 < size <= 16777216              630819            5065558052384
  16777216 < size                        44120533         4342797468343176
 -------------------------------------------------------------------------------
 -                                                                             -
 -                      DBCSR MESSAGE PASSING PERFORMANCE                      -
 -                                                                             -
 -------------------------------------------------------------------------------
 ROUTINE             CALLS      AVE VOLUME [Bytes]
 MP_Group                1
 MP_Bcast               28                     12.
 MP_Allreduce         4230                 325415.
 MP_Alltoall          5355               22325841.
 MP_Wait            200464
 MP_ISend            94336               47211041.
 MP_IRecv            94336               47206344.
 MP_Memory           19172
 -------------------------------------------------------------------------------

 MEMORY| Estimated peak process memory [MiB]                                6402

```

## Key output elements

SCF cycles output looks like
```
 ------------------------------ Linear scaling SCF -----------------------------
 SCF     1  -2019286.740626037  -2019286.740626037  257.119279
 SCF     2  -2022372.558790133     -3085.818164097  283.782252
 SCF     3  -2024495.728609510     -5208.987983474   72.828389
 SCF     4  -2026488.297578288     -7201.556952252  117.651759
 SCF     5  -2029048.217573609     -2559.919995321  229.818745
 SCF     6  -2030970.001095107     -4481.703516819   73.634025
 SCF     7  -2032985.476104558     -6497.178526270  114.296279
 SCF     8  -2033540.457572915      -554.981468357  234.732257
 SCF     9  -2033521.016724448      -535.540619890   74.614194
 SCF    10  -2033612.235695642      -626.759591084  112.382977
 SCF    11  -2033854.276528493      -242.040832851  237.561626
 SCF    12  -2033982.666816269      -370.431120626   66.611201
 SCF    13  -2034003.081946934      -390.846251292   74.652012
 SCF    14  -2034069.478891714       -66.396944780  264.631421
 SCF    15  -2033892.868982235       110.212964699   74.861617
 SCF    16  -2034077.572747632       -74.490800698  108.504719
 SCF    17  -2034120.056681868       -42.483934236  281.836313
 SCF    18  -2034151.889439356       -74.316691724   66.644064
 SCF    19  -2034184.719602434      -107.146854802   74.555730
 SCF    20  -2034205.581420024       -20.861817590  260.906525
 SCF    21  -2034080.981126470       103.738475963   75.491708
 SCF    22  -2034213.969040957       -29.249438523  106.544563
 SCF    23  -2034217.910995587        -3.941954630  318.571562
 SCF    24  -2034219.786293399        -5.817252442   73.672705
 SCF    25  -2034219.952455441        -5.983414485  104.098476
 SCF    26  -2034224.105341386        -4.152885945  295.732344
 SCF    27  -2034225.354788742        -5.402333301   74.348011
 SCF    28  -2034225.361724817        -5.409269375  106.903410
 SCF    29  -2034226.617073632        -1.255348816  312.159937
 SCF    30  -2034226.664080255        -1.302355439   74.832404
 SCF not converged!
 SCF     1  -2034504.114579430      -102.226643770  294.428820
 SCF     2  -2034582.030317720      -180.142382059   59.509095
 SCF     3  -2034669.633673732      -267.745738071   67.661870
 SCF     4  -2034686.756534706       -17.122860974  240.952248
 SCF     5  -2034383.864712474       285.768961258   70.071076
 SCF     6  -2034726.751386791       -57.117713059   97.803430
 SCF     7  -2034732.385823286        -5.634436496  286.235458
 SCF     8  -2034737.205264020       -10.453877229   61.265636
 SCF     9  -2034749.118323179       -22.366936388   69.003241
 SCF    10  -2034748.012737155         1.105586024  305.655801
 SCF    11  -2034751.398910127        -2.280586948   63.232803
 SCF    12  -2034751.430011003        -2.311687824   62.608883
 SCF    13  -2034751.818950539        -0.388939536  333.393301
 SCF    14  -2034752.073886994        -0.643875991   62.769221
 SCF    15  -2034752.205481594        -0.775470591   63.329253
 SCF    16  -2034752.282033039        -0.076551445  351.013058
 SCF    17  -2034752.318892021        -0.113410427   63.425647
 SCF    18  -2034752.322544992        -0.117063398   63.808013
```
