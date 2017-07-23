# Command file to run novaCTF
Algorithm filterProjections
InputProjections corrected_stack_flipped.ali_0
OutputFile filtered_stack.ali_0
TILTFILE angles.tlt
StackOrientation xz


