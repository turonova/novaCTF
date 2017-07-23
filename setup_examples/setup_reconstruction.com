# Command file to run novaCTF
Algorithm 3dctf
InputProjections filtered_stack.ali
OutputFile tomogram.rec
TILTFILE angles.tlt
THICKNESS 140
FULLIMAGE 464 464
SHIFT 0.0 -20.0
PixelSize 0.135
NumberOfInputStacks 10
Use3DCTF 1




