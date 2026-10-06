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

cmake_flags = ["-DCMAKE_BUILD_TYPE=Release", "-DPHX_ACCURACY=FP64"]

if len(sys.argv) > 1 and sys.argv[-1] == "cpu":
    pb_printstr("Building PHOENIX in forced CPU-only mode...")
    cmake_flags.append("-DPHX_FORCE_CPU=ON")
else:
    pb_printstr("Building PHOENIX in default mode...")
    

subprocess.run(["cmake", "-S", ".", "-B", "build"] + cmake_flags, check=True)
subprocess.run(["cmake", "--build", "build", "--", "-j"], check=True)

if os.path.exists("PHOENIX"):
    os.remove("PHOENIX")

if os.path.exists("build/PHOENIX"):
    shutil.move("build/PHOENIX", "./PHOENIX")
    pb_printstr("Build complete! Staged executable at ./PHOENIX")

subprocess.run(["./PHOENIX"], check=True)