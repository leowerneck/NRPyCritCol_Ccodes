import os,sys

def run_string(optimization_string):
    return optimization_string+"./NRPyCritCol 320 2 2 64 0.2 2.293577981651376 0.0"

pulse_center = 0
pulse_width  = 1
Nr           = 200001
R_max        = 70
outputfile   = "outputSFC.txt"

def initial_data_string(eta):
    return "python ScalarField_Initial_Data.py %s %s %lf %lf %d %lf"%(outputfile,eta,pulse_center,pulse_width,Nr,R_max)

def generate_runscript():
    if sys.platform == "linux" or sys.platform == "linux2":

        # Use taskset to improve code's performance
        try:
            import multiprocessing as mp

            # Get only physical cores
            N_cores = int(mp.cpu_count()/2)

            # Write string
            optimization_string = "taskset -c 0"
            for i in range(1,N_cores):
                optimization_string += ","+str(i)
            optimization_string += " "
        except:
            print("Could not import multiprocessing module. Generating runscript without processor affinity.")
            optimization_string += ""

    elif sys.platform == "darwin":

        # Mac OS does not support taskset. Do nothing.
        optimization_string = ""

    else:
        # Windows
        print("Windows detected. runscript.sh will not be generated. Please run the code manually.")
        return

    # Set eta_weak and eta_strong
    eta_weak   = "0.30332394090"
    eta_strong = "0.30332394095"

    # Initialize file with the bash environment
    # This string contains the command to run NRPyCritCol
    runstring   = run_string(optimization_string)+"\n"

    # Initialize the file with the bash environment
    filestring  = "#!/bin/bash\n"

    # Generate weak initial data and make a copy of the ID file
    filestring += initial_data_string(eta_weak)+"\n"
    filestring += "cp %s ID_weak.txt\n"%(outputfile)

    # Run the code; move output to out_weak.dat
    filestring += runstring
    filestring += "mv out_central_values.dat out_weak.dat\n"

    # Generate strong initial data and make a copy of the ID file
    filestring += initial_data_string(eta_strong)+"\n"
    filestring += "cp %s ID_strong.txt\n"%(outputfile)

    # Run the code; move output to out_strong.dat
    filestring += runstring
    filestring += "mv out_central_values.dat out_strong.dat\n"

    # Geerate the plots
    filestring += "python generate_plot.py"

    with open("runscript.sh","w") as file:
        # Write string to file
        file.write(filestring)

    os.system("chmod 755 runscript.sh")

if __name__ == '__main__':
    generate_runscript()
