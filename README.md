# R code for fitting penalized multistate additive hazards regression models

## Overview 

This repository contains an R-based implementation of the methodology presented in the following submitted paper:

*Penalized Variable Selection in Additive Hazards Regression for Multi-State Time-to-Event Data: Application to COVID-19*

Note: The implementation is based on and extends the original code by Dr. Jinchi Lv, available at the Marshall School of Business, University of Southern California: [Lv's Software Page](https://faculty.marshall.usc.edu/jinchi-lv/publications/software/)


## Table of Contents

* [Included Files](#includedfiles)
* [Compilation Notes](#compilationnotes)
* [Usage](#usage)
* [Troubleshooting](#troubleshooting)
* [Final Remarks](#finalremarks)
* [Reference](#reference)


## Included Files

* `ahfun.c`              – C source code  
* `ahfun.R`              – R wrapper  
* `Simulation.example.R` – Example R script demonstrating usage  
* `ahfun32.dll`          – DLL (Dynamic-Link Library) compiled from `ahfun.c` for Windows 32-bit
* `ahfun64.dll`          – DLL (Dynamic-Link Library) compiled from `ahfun.c` for Windows 64-bit 
* `README`               – Project documentation (this file)

Make sure that the file `ahfun.R` and one of the following files (depending on your operating system) are available in R's working directory:

* `ahfun32.dll`
* `ahfun64.dll`
* `ahfun.so`          – Shared object (SO) compiled from `ahfun.c` for Linux
 
 Note: On Linux, you must manually compile `ahfun.c` to create `ahfun.so`. Please refer to the [Compilation Notes](#compilationnotes) section for instructions.


## Compilation Notes

* **Windows**: You can directly use the appropriate precompiled DLL (`ahfun32.dll` or `ahfun64.dll`) based on your system architecture. Next, run the `Simulation.example.R` script.

* **Linux**: Before executing `Simulation.example.R`, you must first compile the C source code by entering the following command in the terminal:  
  
  R CMD SHLIB ahfun.c

This command will generate a shared object file named `ahfun.so`.


Note: If you are working on Windows and want to compile the source code yourself, make sure Rtools (https://cran.r-project.org/bin/windows/Rtools) is installed.


## Usage

The `Simulation.example.R` file contains the main implementation code, along with inline comments and section headers that guide the user through installing necessary packages, setting up the data, running the penalized multistate additive hazards regression model, and saving the results.


## Troubleshooting

* DLL File Errors (Windows): If you encounter an error related to a DLL file while running `Simulation.example.R` on Windows, make sure you are using the correct version (`ahfun32.dll` or `ahfun64.dll`) that matches your system's architecture (Windows 32-bit or Windows 64-bit).

* Compilation Issues (Linux): If you experience problems compiling the C source code on Linux, make sure that essential development tools (such as gcc and make) are properly installed on your system.


## Final Remarks

Thank you for using this repository. If you encounter any issues or have questions, feel free to open an issue on GitHub, and we will be happy to assist you. Contributions are always welcome! We hope this implementation proves useful for your research.


## Reference

* Tapak, L., Bandyopadhyay, D., Hamidi, O., & Arayeshgari, M. (2025+). *Penalized variable selection in additive hazards regression for multi-state time-to-event data: Application to COVID-19*, (submitted).

* Lin, W., & Lv, J. (2013). *High-dimensional sparse additive hazards regression*. Journal of the American Statistical Association, 108(501), 247–264.
