f2py -c -m chapman chapman.for

# Or, for only a few functions,
f2py -c chapman.for only: atm8_chapman_arr atm8_chapman atm8_chap_num : -m chapman

#Or, for stuff that should let you use 'intrinsic ieee_value' but doesn't actually,
f2py -c chapman.for --opt="-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans" only: atm8_chapman_arr atm8_chapman atm8_chap_num : -m chapman
