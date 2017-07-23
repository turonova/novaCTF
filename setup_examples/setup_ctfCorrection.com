# Command file to run novaCTF
Algorithm ctfCorrection
InputProjections input_stack.st
DefocusFile defocus_file.txt_0
OutputFile corrected_stack.st_0
TILTFILE angles.tlt
CorrectionType phaseflip 
DefocusFileFormat ctffind4
PixelSize 0.135
AmplitudeContrast 0.07
Cs 2.7
Volt 300
CorrectAstigmatism 1


