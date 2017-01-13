#__GPS2TEC__

This prohram is for calculate total electron content(TEC) from gps signal



#Directory and file discription:

__bin/__  : some binary files (different version of crx2rnx)

__gps2tec_modules/__  : subroutines



__gps2tec.inp__ :  input file 

__gps2tec.py__  : main program 

__marker.crd__  : gps station xyz coordinate list 

#Input file discription:
__case_type__  : local  -- Calculate GPSTEC from o file in current directory
                 IGS    -- Download GPS source file from IGS(International GNSS Service) and calculate TEC
                 IGSRT  -- Download real-time GPS source file from IGS and calculate TEC
