---------
OVERVIEW:
---------

This directory contains a program that calculates parameters for the F82-Tint metal reflectivity model for several common and interesting metals. The calculated parameters can be directly plugged into the Adobe Standard Material (ASM) or OpenPBR material models, which both use the F82-Tint model.

The program is intended for experimentation, validation, and comparison. The source code is rough and pragmatic and is not intended to be used in production. Since the calculations already complete very quickly, the code does not need to be optimized. All code is kept in a single file for convenience.

---------------------
BUILDING AND RUNNING:
---------------------

To build and run the program, navigate to this directory and run this command:

./BuildAndRun.sh

Input data is automatically read from the various data files at runtime. To add more metals or upgrade the data for any metals, add more files to the ComplexRefractiveIndexes directory.

The program will output the results into the OpenPBRBaseColorAndSpecularColorForRealMetals.txt file in this directory.

------------------
TECHNICAL DETAILS:
------------------

The F82-Tint material model is an extension of the F82 model. It was introduced and described in the ASM technical documentation here:

https://helpx.adobe.com/substance-3d-general/adobe-standard-material/asm-technical-documentation.html

The parameters output by the program can be used in ASM and OpenPBR to set the following parameters:

ASM:

- base_color
- specular_edge_color

OpenPBR:

- base_color
- specular_color
