import sys
import os
import numpy as np


XHI_QSO_Fan2006_z = np.array([5.00, 5.25, 5.44, 5.65])
XHI_QSO_Fan2006_XHI = np.log10(np.array([5.49e-5, 6.75e-5, 6.50e-5, 8.56e-5]))
XHI_QSO_Fan2006_XHImax = np.log10(np.array([7.00e-5, 8.91e-5, 9.21e-5, 1.35e-4]))
XHI_QSO_Fan2006_XHImin = np.log10(np.array([4.02e-5, 4.45e-5, 4.14e-5, 4.75e-5]))

XHI_QSO_Fan2006_lowlim_z = np.array([5.84, 6.10])
XHI_QSO_Fan2006_lowlim_XHI = np.log10(np.array([1.25e-4, 4.36e-4]))
XHI_QSO_Fan2006_lowlim_XHImax = np.log10(np.array([3.e-4, 1.e-3]))
XHI_QSO_Fan2006_lowlim_XHImin = np.log10(np.array([7.699e-5, 1.338e-4]))


XHI_GRB_Totani2014_z = np.array([5.92])
XHI_GRB_Totani2014_XHI = np.log10(np.array([0.298]))
XHI_GRB_Totani2014_XHImax = np.log10(np.array([0.501]))
XHI_GRB_Totani2014_XHImin = np.log10(np.array([0.102]))

XHI_GRB_Totani2014_uplim_z = np.array([6.3])
XHI_GRB_Totani2014_uplim_XHI = np.log10(np.array([0.166]))


XHI_LyaLF_Konno2018_z = np.array([6.6, 7.3])
XHI_LyaLF_Konno2018_XHI = np.log10(np.array([0.302, 0.555]))
XHI_LyaLF_Konno2018_XHImax = np.log10(np.array([0.501, 0.802]))
XHI_LyaLF_Konno2018_XHImin = np.log10(np.array([0.102, 0.302]))

XHI_LyaLF_Malhotra2004_z = np.array([6.5])
XHI_LyaLF_Malhotra2004_uplim_XHI = np.log10(np.array([0.30]))

XHI_LyaLF_Ota2010_z = np.array([7.0])
XHI_LyaLF_Ota2010_uplim_XHI = np.log10(np.array([0.63]))

XHI_LyaLF_Ouchi2010_z = np.array([6.7])
XHI_LyaLF_Ouchi2010_uplim_XHI = np.log10(np.array([0.20]))

XHI_LyaLF_Kashikawa2011_z = np.array([6.55])
XHI_LyaLF_Kashikawa2011_uplim_XHI = np.log10(np.array([0.40]))


XHI_LAEACF_Ouchi2010_z = np.array([6.76])
XHI_LAEACF_Ouchi2010_uplim_XHI = np.log10(np.array([0.501]))


XHI_LAEfrac_z = np.array([7.12])
XHI_LAEfrac_lowlim_XHI = np.log10(np.array([0.501]))


tau = 0.0544
tau_uplim = 0.0544 + 0.0070
tau_lowlim = 0.0544 - 0.0081