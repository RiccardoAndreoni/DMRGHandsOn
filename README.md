# DMRGHandsOn
This repository contains the code to run a simple DMRG simulation for a Heisenberg model on a chain.

## How to compile
Use make. In order to compile all executables do:  
```bash
mkdir build
cd build
cmake ..
make
```

In order to delete all compiled files do:  
```bash
make clean
```

## Structure
```bash
.
├── include                 # Headers
├── modules                 # Routines
├── main                    # Source files for the main programs
├── test                    # Source files for the test programs
├── CMakelists.txt
├── .gitignore
└── README.md
```