# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 17:02:47 2023

@author: jackd
"""

## SD errors;

sd_sample01_reserrs = pd.read_csv('./SD_sample01_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
sd_sample02_reserrs = pd.read_csv('./SD_sample02_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
sd_sample03_reserrs = pd.read_csv('./SD_sample03_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
sd_sample04_reserrs = pd.read_csv('./SD_sample04_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)

sd_sample01errs = np.array(sd_sample01_reserrs.iloc[0])
sd_sample02errs = np.array(sd_sample02_reserrs.iloc[0])
sd_sample03errs = np.array(sd_sample03_reserrs.iloc[0])
sd_sample04errs = np.array(sd_sample04_reserrs.iloc[0])

sd_diff2banded10_reserrs = pd.read_csv('./SD_diff2banded10_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
sd_diff2banded25_reserrs = pd.read_csv('./SD_diff2banded25_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
sd_diff2banded50_reserrs = pd.read_csv('./SD_diff2banded50_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)

sd_diff2banded10errs = np.array(sd_diff2banded10_reserrs.iloc[0])
sd_diff2banded25errs = np.array(sd_diff2banded25_reserrs.iloc[0])
sd_diff2banded50errs = np.array(sd_diff2banded50_reserrs.iloc[0])

sd_iters01 = len(sd_sample01errs)
sd_iters02 = len(sd_sample02errs)
sd_iters01 = np.array( range(len(sd_sample01errs)) )
sd_iters02 = np.array( range(len(sd_sample02errs)) )
sd_iters03 = np.array( range(len(sd_sample03errs)) )
sd_iters04 = np.array( range(len(sd_sample04errs)) )
sd_iters10 = np.array( range(len(sd_diff2banded10errs)) )
sd_iters25 = np.array( range(len(sd_diff2banded25errs)) )
sd_iters50 = np.array( range(len(sd_diff2banded50errs)) )

plt.xlabel("number of iterations"); \
plt.ylabel("log10 residual errors"); \
plt.title("SD errors"); \
plt.plot(sd_iters01, np.log10(sd_sample01errs), sd_iters02, np.log10(sd_sample02errs), sd_iters03, np.log10(sd_sample03errs), sd_iters04, np.log10(sd_sample04errs));

plt.xlabel("number of iterations"); \
plt.ylabel("log10 residual errors"); \
plt.title("SD errors"); \
plt.plot( sd_iters10, np.log10(sd_diff2banded10errs), sd_iters25, np.log10(sd_diff2banded25errs), sd_iters50, np.log10(sd_diff2banded50errs) );

## CGD errors;

cgd_sample01_reserrs = pd.read_csv('./CGD_sample01_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
cgd_sample02_reserrs = pd.read_csv('./CGD_sample02_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
cgd_sample03_reserrs = pd.read_csv('./CGD_sample03_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
cgd_sample04_reserrs = pd.read_csv('./CGD_sample04_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)

cgd_sample01errs = np.array(cgd_sample01_reserrs.iloc[0])
cgd_sample02errs = np.array(cgd_sample02_reserrs.iloc[0])
cgd_sample03errs = np.array(cgd_sample03_reserrs.iloc[0])
cgd_sample04errs = np.array(cgd_sample04_reserrs.iloc[0])

cgd_diff2banded10_reserrs = pd.read_csv('./CGD_diff2banded10_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
cgd_diff2banded25_reserrs = pd.read_csv('./CGD_diff2banded25_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)
cgd_diff2banded50_reserrs = pd.read_csv('./CGD_diff2banded50_results.txt', skiprows=4, nrows=1, delim_whitespace=True, header=None)

cgd_diff2banded10errs = np.array(cgd_diff2banded10_reserrs.iloc[0])
cgd_diff2banded25errs = np.array(cgd_diff2banded25_reserrs.iloc[0])
cgd_diff2banded50errs = np.array(cgd_diff2banded50_reserrs.iloc[0])

cgd_iters01 = len(cgd_sample01errs)
cgd_iters02 = len(cgd_sample02errs)
cgd_iters01 = np.array( range(len(cgd_sample01errs)) )
cgd_iters02 = np.array( range(len(cgd_sample02errs)) )
cgd_iters03 = np.array( range(len(cgd_sample03errs)) )
cgd_iters04 = np.array( range(len(cgd_sample04errs)) )
cgd_iters10 = np.array( range(len(cgd_diff2banded10errs)) )
cgd_iters25 = np.array( range(len(cgd_diff2banded25errs)) )
cgd_iters50 = np.array( range(len(cgd_diff2banded50errs)) )

plt.xlabel("number of iterations"); \
plt.ylabel("log10 residual errors"); \
plt.title("CGD errors"); \
plt.plot(cgd_iters01, np.log10(cgd_sample01errs), cgd_iters02, np.log10(cgd_sample02errs), cgd_iters03, np.log10(cgd_sample03errs), cgd_iters04, np.log10(cgd_sample04errs));

plt.xlabel("number of iterations"); \
plt.ylabel("log10 residual errors"); \
plt.title("CGD errors"); \
plt.plot( cgd_iters10, np.log10(cgd_diff2banded10errs), cgd_iters25, np.log10(cgd_diff2banded25errs), cgd_iters50, np.log10(cgd_diff2banded50errs) );

## flatlevel function;

def flatlevel(errs_list):
"""
flatten zero values to 10e-7;
"""
    for i,x in enumerate(errs_list):
        errs_list[i] = 10e-7 if x==0.0 else x;
## ~ fin flatlevel ##