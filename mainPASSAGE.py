import wisp_analysis as wisp
import os

# Open two ds9 windows:
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_DIRECT &')
os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_spec2D &')

CODE_DIR = "/Users/knedkova/Work/2024PASSAGE/wisp_analysis"
OUTPUT_DIR = "/Users/knedkova/Work/2024PASSAGE/output"
path_to_wisp_data = "/Users/knedkova/Work/2024PASSAGE/data/"

if __name__ == "__main__":
    # modify the spec2D directory to stamps...
    os.chdir(OUTPUT_DIR)
    # wisp.loop_field_cwt()
    wisp.measure_z_interactive(path_to_wisp_data=path_to_wisp_data, path_to_code=CODE_DIR)

