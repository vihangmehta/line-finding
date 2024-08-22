import wisp_analysis as wisp
import os

# Open two ds9 windows:
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_DIRECT &')
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_spec2D &')

CODE_DIR = "/Users/knedkova/Work/2024PASSAGE/wisp_analysis"
OUTPUT_DIR = "/Users/knedkova/Work/2024PASSAGE/output"
path_to_wisp_data = "/Users/knedkova/Work/2024PASSAGE/data/"

#############################################################
####### You should not need to change anything below. #######
#############################################################

parno = input('\033[94m' + "Enter the number of the parallel field you want to analyze.\n> " + '\033[0m')
# pull out only the number, in case the user entered e.g. "Par1"
while True:
    try:
        parno = int(re.findall(r'\d+', str(parno))[0])
    except: 
        parno = input('\033[94m' + "A parallel field number is required. Enter the number of the parallel field you want to analyze.\n> " + '\033[0m')
        continue
    else:
        break

if __name__ == "__main__":
    # move to the directory where you want your outputs.
    os.chdir(OUTPUT_DIR)

    # check if region files exist. If not, run code to create necessary region files
    regionfiles = glob.glob(path_to_wisp_data + "Par" + str(parno) + "/DATA/*.reg")
    if len(regionfiles) == 0:
        print('\033[94m' + "No region files found, creating those for you now."  + '\033[0m')
        create_regions(parno = parno, path_to_wisp_data = path_to_wisp_data)

    # check if line list exists. If not, run code to create the linelist
    print(OUTPUT_DIR + "/linelist/Par"+str(parno)+"lines.dat")
    linelist = glob.glob(OUTPUT_DIR + "/linelist/Par"+str(parno)+"lines.dat")
    if len(linelist) == 0:
        print('\033[94m' + "No line list file found, creating the line list for you now." + '\033[0m')
        wisp.loop_field_cwt(path_to_wisp_data=path_to_wisp_data, path_to_code=CODE_DIR, parno=parno)
    
    # run the measure_z_interactive codes
    wisp.measure_z_interactive(path_to_wisp_data=path_to_wisp_data, path_to_code=CODE_DIR, parno=parno)


