import numpy as np
from grid import *

def compute_cross_corr(field1, field2, boxsize, outputfile, str1, str2):

    field1_norm = field1/np.mean(field1, dtype=np.float64) - 1.
    field2_norm = field2/np.mean(field2, dtype=np.float64) - 1.

    modesField1 = ifftn(field1_norm)
    modesField2 = ifftn(field2_norm)

    kmid_bins, ps_field1, ps_field1_err = modes_to_pspec(modesField1, boxsize=boxsize)
    kmid_bins, ps_field2, ps_field2_err = modes_to_pspec(modesField2, boxsize=boxsize)
    kmid_bins, crossps, crossps_err = two_modes_to_pspec(modesField1, modesField2, boxsize=boxsize)

    crosscorreff = crossps / (ps_field1 * ps_field2)**0.5

    sigma_b = ps_field1
    sigma_c = ps_field2
    sigma_a2 = 0.5 * (crossps**2 + sigma_b * sigma_c)
    sigma_ab = crossps * ps_field1
    sigma_ac = crossps * ps_field2
    sigma_bc = ps_field1 * ps_field2

    tmp = sigma_a2 / crossps**2 + (0.5 * sigma_b / ps_field1)**2 + (0.5 * sigma_c / ps_field2)**2 - sigma_ab / (crossps * ps_field1) - sigma_ac / (crossps * ps_field2) + 0.5 * sigma_bc / (ps_field1 * ps_field2)

    crosscorreff_err = crosscorreff * tmp**0.5

    bias = ps_field1 / ps_field2

    np.savetxt(outputfile + 'ps_' + str1 + '.dat', np.c_[kmid_bins, ps_field1, ps_field1_err])
    np.savetxt(outputfile + 'ps_' + str2 + '.dat', np.c_[kmid_bins, ps_field2, ps_field2_err])
    np.savetxt(outputfile + 'crossps_' + str1  + '_' + str2 + '.dat', np.c_[kmid_bins, crossps, crossps_err])
    np.savetxt(outputfile + 'crosscorreff_' + str1  + '_' + str2 + '.dat', np.c_[kmid_bins, crosscorreff, crosscorreff_err])
    np.savetxt(outputfile + 'bias_' + str1 + '_' +str2 + '.dat', np.c_[kmid_bins, bias])

