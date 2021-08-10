#ifdef __cplusplus
extern "C" {
#endif
struct NNP;
typedef struct NNP nnp;

nnp *create_nnp(char *model);

void delete_nnp(nnp *n);

void compute_nnp(nnp *n, int *vecsize, double *dener, double *dforce,
                 double *dvirial, double *datom_ener, double *datom_virial,
                 double *dcoord_, int *datype_, double *dbox);

#ifdef __cplusplus
}
#endif
