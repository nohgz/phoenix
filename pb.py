#
# pb.py (PyBuild.py)
# ---
# A convenience script that runs the annoying CMake stuff for you!
#

import os
import shutil
import subprocess
import sys

def pb_printstr(str):
    print("[\033[1m\033[38;5;44mPyBuild\033[0m] " + str)

shutil.rmtree("build", ignore_errors=True)

args = sys.argv[1:]
ACCURACY = ("FP32", "FP64")
accuracy = next((a.upper() for a in args if a.upper() in ACCURACY), "FP32")
cmake_flags = ["-DCMAKE_BUILD_TYPE=Release", "-DPHX_ACCURACY=" + accuracy]

if "cpu" in args:
    pb_printstr("Building PHOENIX in forced CPU-only mode...")
    cmake_flags.append("-DPHX_FORCE_CPU=ON")
else:
    pb_printstr("Building PHOENIX in default mode...")

if "phi" in args:
    pb_printstr("Enabling azimuthal camera!") # this is kind of dumb, and im unsure what kind of gains it warrants. why did i do this
    cmake_flags.append("-DPHX_CAMERA_PHI=ON")

subprocess.run(["cmake", "-S", ".", "-B", "build"] + cmake_flags, check=True)
subprocess.run(["cmake", "--build", "build", "--", "-j"], check=True)

if os.path.exists("PHOENIX"):
    os.remove("PHOENIX")

if os.path.exists("build/PHOENIX"):
    shutil.move("build/PHOENIX", "./PHOENIX")
    pb_printstr("Build complete! Staged executable at ./PHOENIX")

if "run" in args:
    subprocess.run(["./PHOENIX"], check=True)