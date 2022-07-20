# cgem2amrvac

This version so far only runs on intel chipsets

Required to run this pipeline:
 - f2py(numpy)
 - sdf
 - fishpack
 - PDFI

files included and their purpose
 - pipeline.ipynb           - contains the python script required to download CGEM_input files from JSOC given a sharp number and a time range
 - save_pdfi_to_sdf.f95     - contains the relevant calls from the test_wrapper function in the PDFI software

One must first compile the save_pdfi_to_sdf.f95 file using f2py following:

f2py -c --opt="-fallow-argument-mismatch -fcheck=all -fbounds-check" save_pdfi_to_sdf.f95 -m save_pdfi_to_sdf -L/usr/local/lib/ -lpdfi_ss -lsdf -lfishpack

or

f2py -c --opt="-fcheck=all -fbounds-check" save_pdfi_to_sdf.f95 -m save_pdfi_to_sdf -L/users/cpa/jackj/codes/dd_lib/lib -lpdfi_ss -lsdf -lfishpack

This subroutine may then be called from within python - see pipeline.ipynb


Sharp ARs known to have a CGEM counterpart, and also successfully passed from python to fortran, imported into amrvac and simulated:
 - 377
 - 401
