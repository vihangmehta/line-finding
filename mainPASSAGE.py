import wisp_analysis as wisp
import os

# Open two ds9 windows:
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_DIRECT &')
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_spec2D &')

CODE_DIR = "/Users/knedkova/Work/2024PASSAGE/wisp_analysis"
OUTPUT_DIR = "/Users/knedkova/Work/2024PASSAGE/output"
path_to_wisp_data = "/Users/knedkova/Work/2024PASSAGE/data/"

# The field number:
parno = 28

# To run the code, first run wisp.loop_field_cwt() to create the line list for thr field you're working on
# Once that is done, that line can be commented out, unless the line list needs to be recreated (e.g. if you're changing fields).
# Next, run wisp.measure_z_interactive() to fit the spectra.

if __name__ == "__main__":
    # move to the directory where you want your outputs.
    os.chdir(OUTPUT_DIR)
    # wisp.loop_field_cwt(path_to_wisp_data=path_to_wisp_data, path_to_code=CODE_DIR, parno=parno)
    wisp.measure_z_interactive(path_to_wisp_data=path_to_wisp_data, path_to_code=CODE_DIR, parno=parno)

