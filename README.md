Control Systems toolbox for Julialang

Implementation plan:
=====

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

    [x] State Space

        [x] Constructors
        [x] Printing
        [x] Addition
        [x] Subtraction
        [x] Multiplication
        [x] Division
            -NOTE: Currently division by SS objects isn't supported. However, this is the current state with python-control as well. Low key issue. Will examine Octave's control lib at somepoint to implement this feature.

    [ ] tf2ss
    [x] ss2tf
            -NOTE: This depends on Slicot.jl. Haven't figured out packaging yet, so this isn't robust at all. TODO.

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
