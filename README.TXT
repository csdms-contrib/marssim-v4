The package includes 6 directories:
    accessory_programs:  A variety of programs to create input files, particularly pseudo-fractal surfaces, to
        extract a specified matrix file from a concatinated output file, such as OUTELEV.DAT.  Also a program
        tosurfer.F to reformat an output file to a form compatible with the commercial analysis program SURFER.
    documentation:  The documetation marssim_program_v4.doc file and the GUI interface with its subdirectories
    execution_directory_files:  the parameter files, example files for a simple simulation, clearerode.sh (linux)
        and clearerode.(windows) erases any lingering files from previous runs. runname.sh (linux)
        runname.bat (windows) renames output files to have a six letter/number prefix (from the command line)
        and creates a zipped archive. ERODE.STOP.NO, ERODE.STOP, and ERODE.STOP.YES for controlling program execution.
        The parameter files, if executed with a compiled version of MARSSIM, will eventually progress toward a steady-
        state fluvial landform using a combination of mass wasting and fluvial downcutting with a constant rate of lowering
        of the lower boundary.
    parameter_files:  Another copy of the default parameter files
    real_craters:  A directory containing the files necessary to simulate new craters using martian fresh craters
    source_files:  The source files for MARSSIM plus example compilation and files for linux using intel fortran.  
        This uses the –fast option which requires static clib libraries. If these are not installed substitute
        –O3 for -fast.  For debugging substitute –g for –fast.
