#
# This script works only with Python 3
#
import os
import re

import create_initial_conditions
import Test

MaxParticlesInSequentialUpscalingStudies = 200


def step1(compiler, default_flags):
  test = Test.Test( compiler, default_flags, "step-1" )
  test.cleanall()

  arguments = "-O3 --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  _,arguments = create_initial_conditions.create_grid_setup( 2,1,1, 1, 1, "no-noise" , 0.1, 0.1, 1 )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()

  search_pattern = "([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+)"
  result = test.search_for_pattern_in_output(search_pattern)
  if result!="":
    print( "Got " + result + " as last output which I interprete to be a particle position ... ok (though that does not mean that the data is correct; that's something I don't validate here)" )
  else:
    print( "Last line of output should still be the one I used in the template ... failed" )
    exit()

def step2(compiler, default_flags):
  test = Test.Test( compiler, default_flags, "step-2" )

  arguments = "-O3 --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  _,arguments = create_initial_conditions.create_grid_setup( 2,1,1, 1, 1, "no_noise" , 0.1, 0.1, 1 )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok (but this does not mean that the outcome is correct)" )
  else:
    print( "Run failed: " + test.last_output )
    exit()

def step3(compiler, default_flags):
  test = Test.Test( compiler, default_flags, "step-3" )
  test.cleanall()

  arguments = "-O0 -fno-tree-vectorize"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [13,12,12]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  no_vec_time = test.runtime
  test.cleanall()

  arguments = "-O3 --std=c++0x -fno-math-errno"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [13,12,12]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  with_vec_time = test.runtime

  if no_vec_time<=with_vec_time:
    print( "Code is slower with vectorisation, so you might want to tune it ... check" )
    exit()
  else:
    print( "Code is already faster by a factor of " + str(no_vec_time/with_vec_time) + " through vectorisation but you might want to tune it further ... ok" )


def step4(compiler, default_flags):
  test = Test.Test( compiler, default_flags, "step-4" )

  arguments = "-O3 --std=c++0x -fno-math-errno"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )
    exit()

  particles_counts = [13,12,12]
  _,arguments = create_initial_conditions.create_grid_setup( particles_counts[0],particles_counts[1],particles_counts[2], 1, 1, "random" , 0.1, 0.0001, 10 )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )

  print( "Test for one core" )
  result = test.run( arguments, {"OMP_NUM_THREADS": "1"} )
  if result==1:
    print( "Run terminated successfully after " + str(test.runtime) + "s ... ok" )
  else:
    print( "Run failed: " + test.last_output )
    exit()
  serial_runtime = test.runtime

  for p in range(2,28,2):
    print( "Test for " + str(p) + " cores" )
    result = test.run( arguments, {"OMP_NUM_THREADS": str(p)} )
    speedup = serial_runtime/test.runtime
    if result==1 and speedup>p*0.8:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. with speedup of " + str(speedup) + " ... ok  (but you might want to tune it further)" )
    elif result==1:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. without expected speedup ... check runtimes" )
    else:
      print( "Run failed: " + test.last_output )
      exit()


def run_all_steps(compiler, default_flags):
  try:
    print("Checking step 1")
    step1(compiler, default_flags)
    print("Checking step 2")
    step2(compiler, default_flags)
    print("Checking step 3")
    step3(compiler, default_flags)
    print("Checking step 4")
    step4(compiler, default_flags)

  except BaseException as err:
    print(f"Unexpected {err=}, {type(err)=}")
    raise


if __name__=="__main__":
  print("""
 This is a small Python script to check that the format of a submission is
 valid. It does not check for correctness, although it gives some clues on
 performance. If you use the script on Hamilton, which is the machine
 we will use to mark the submission, you need to load appropriate modules:

 module load python/3.10.8 gcc/12.2

 We recommended that you use a compute node rather than the login nodes. You
 can do this with the command:

 salloc -N 1 -p test.q python3 validate.py

 Disclaimer: This is only a sanity check, and we'll run further tests on
 correctness and scalability to mark the coursework. But the script ensures
 that your submission is formally correct and performs some basic checks.

""")
  run_all_steps("g++", "-fopenmp -march=native")
