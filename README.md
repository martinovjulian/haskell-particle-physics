## Description 
This project was derived from https://www.cs.columbia.edu/~sedwards/classes/2021/4995-fall/reports/FluidDyn.pdf. The original code had a lot of issues/compiler warnings and didn't follow idiomatic Haskell particularly well. It also did not function on Windows at all. I have successfully ported the code which now works on Windows 11(and should work on any OS).

## Build & Run Instructions
- You MUST have GHC and Stack installed
- Clone the repo
- Run "stack build"
    - This should install all necessary dependencies and build an executable
- Run "stack exec simulation-exe" to run the simulation executable

This is a WIP.