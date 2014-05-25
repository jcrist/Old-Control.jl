# Status Page

## Implementation plan:

Below follows a list of things to do. As they get done, they will be checked off.

Types
---
    [x] Transfer Function

        [x] Constructors
        [x] Printing
        [x] Addition
        [x] Subtraction
        [x] Multiplication
        [x] Division
            -NOTE: Currently division by MIMO TransferFunctions isn't supported.
            However, this is the current state with python-control as well. 
            With addition of LU-Decomposition for numeric types without pivoting
            (in the works), this will be supported.

    [x] State Space

        [x] Constructors
        [x] Printing
        [x] Addition
        [x] Subtraction
        [x] Multiplication
        [x] Division

    [ ] tf2ss
    [ ] ss2tf

Analysis
---
    [ ] pole
    [ ] zero
    [ ] zpkdata
    [ ] tfdata
    [ ] ctrb
    [ ] obsv

Plots
---
    [ ] rlocus
    [ ] pzmap
    [ ] bode
    [ ] step
    [ ] nyquist
    [ ] impuls
    [ ] lsim

Control
---
    [ ] pid
    [ ] feedback

Design
---
    [ ] place
    [ ] acker
    [ ] lqr
    [ ] lqi

Digital
---
    [ ] c2d
