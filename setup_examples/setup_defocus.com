# Command file to run novaCTF
Algorithm defocus
InputProjections input_stack.st 
FULLIMAGE 464 464
THICKNESS 140
TILTFILE angles.tlt
SHIFT 0.0 0.0
CorrectionType phaseflip
DefocusFileFormat ctffind4
DefocusFile defocus_file.txt
PixelSize 0.135
DefocusStep 15
CorrectAstigmatism 1
DefocusShiftFile file_with_additional_defocus.txt

