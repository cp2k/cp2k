import numpy as np
from dftd3.interface import RationalDampingParam, DispersionModel

num = np.array([6, 1, 1, 1, 1])
xyz = np.array(  # coordinates in Bohr
    [
        [ 0.0000000, -0.0000000,  0.0000000],
        [-1.1922080,  1.1922080,  1.1922080],
        [ 1.1922080, -1.1922080,  1.1922080],
        [-1.1922080, -1.1922080, -1.1922080],
        [ 1.1922080,  1.1922080, -1.1922080],
    ]
)
method = "PBE0"

model = DispersionModel(num, xyz)
res = model.get_dispersion(RationalDampingParam(method=method), grad=False)
print(f"Dispersion energy for {method}-D3(BJ) is {res['energy']:13.10f} Hartree")
