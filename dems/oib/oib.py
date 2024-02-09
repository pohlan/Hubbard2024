import os
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description='PDAL DEM generation pipeline wrapper.')
parser.add_argument("input", help="Dir of .las files to process")
parser.add_argument("output", help="Dir to save .tif output to")
parser.add_argument("-r", "--resolution", help="Pixel resolution of output DEM (Default = 10)", default=10, type=float)
args = parser.parse_args()

files = glob.glob("%s/*.las" % args.input)

fd = open("jobs.txt", "w")

for file in files:
	bn = os.path.basename(file)
	out = args.output + "/" + bn.replace(".las", "_%.0fm.tif" % args.resolution)
	fd.write("pdal pipeline oib.json --readers.las.filename=%s --writers.gdal.resolution=%f --writers.gdal.filename=%s\n" % (file, args.resolution, out))

fd.close()
