# Fast-Fourier-Transform ⚡

In this project I implement the Fast Fourier Transform algorithm using C++. Moreover I also using other optimization techniques such as strength reduction to produce fast and efficicent audio convolution.

## How to run my program :tada:

There are two different folders in my repository named “base_program” and “FFT_program”, the code for the base program along with its output and gprof file are located in the base_program folder. Similarly, the code for the FFT program along with all the different outputs for each optimization and all the gprofs for each optimization are located in the “FFT_program” folder. I have made all optimizations on the FFT.cpp file itself so in order to see unoptimized and optimize progression you will have to look at the commit history. For both (base.cpp and FFT.cpp) programs compilation is identical. Here is an example:

    g++ base.cpp
    ./a.out guitar_dry.wav big_hall_IR_mono.wav output.wav
    
    g++ FFT.cpp
    ./a.out guitar_dry.wav big_hall_IR_mono.wav output.wav
