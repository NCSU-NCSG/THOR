
# Copyright (c) 2005-2010, 2012-2013, Andrew Hang Chen and contributors,
# All rights reserved.
# Licensed under the 3-clause BSD license.

module RakeBase
  puts "RUBY_PLATFORM=" + RUBY_PLATFORM if $show_info

  if RUBY_PLATFORM =~ /(darwin|linux)/i
    # Intel FORTRAN compiler tested on Linux
    $compiler = 'ifort'
    $option = "-check all -warn all -fpp"
    $ext_obj = "o"
    $dosish_path = false
    $gcov = false
		$prof_genx = "-prof-genx"
    $mpiexec = nil
  else
    # Intel FORTRAN on Windows
    $compiler = 'ifort'
    $option = "/check:all /warn:all /fpp"
    $ext_obj = "obj"
    $dosish_path = true
    $gcov = false
		$prof_genx = "/Qprof-genx"
    $mpiexec = nil
  end


  # GCC FORTRAN compiler tested on MacOs (10.6.8 Snow Leopard) and Windows Vista + cygwin
  #$compiler = "gfortran"
  #$option = "-Wall -Wextra -pedantic -fbounds-check " +
  #          "-Wuninitialized -O -g -Wno-unused-parameter -cpp"
  #$ext_obj = "o"
  #$dosish_path = false
  # # With " -std=f95",
  # # subroutines whose name is longer than 31 characters cause error.


  # #G95 FORTRAN compiler tested on Linux and Windows Vista + cygwin
  #$compiler = "g95"
  #$ext_obj = "o"
  #$dosish_path = false
  #$option = "-Wall -Wobsolescent -Wunused-module-vars " +
  #  "-Wunused-internal-procs -Wunused-parameter -Wunused-types " +
  #  "-Wmissing-intent -Wimplicit-interface -pedantic -fbounds-check -Wuninitialized"


  # #FTN95 Fortran compiler
  #  $compiler = "ftn95"
  #  $linker = "slink"
  #  $option = "/CFPP /DEFINE FTN95 1 /SILENT "
  #  $linker_option = ""
  #  $option_obj = "/binary"
  #  $ext_obj = "obj"
  #  $option_exe = "-out:"
  #  $dosish_path = true

  #---------------------------------------------------
  #`where ...` works on windows vista, 7 and 8. Not works on Windows XP.
  test_where = `where where 2>&1`
  if $?.to_i == 0
    where_or_which = "where"
  else
    where_or_which = "which"
  end

  #----- if absent try gfortran
  result = `#{where_or_which} #{$compiler} 2>&1`
  if $?.to_i != 0
    puts "Fortran compiler " + $compiler + " not exists. Trying gfortran."
    $compiler = "gfortran"
    $option = "-Wall -Wextra -pedantic -fbounds-check " +
              "-Wuninitialized -O -g -Wno-unused-parameter -cpp "
    $ext_obj = "o"
    $dosish_path = false
    $gcov = "-coverage"
		$prof_genx = false
  end

  # ----- if absent try FTN95
  result = `#{where_or_which} #{$compiler} 2>&1`
  if $?.to_i != 0
    puts "Fortran compiler " + $compiler + " not exists. Trying ftn95."
    $compiler = "ftn95"
    result = `where #{$compiler} 2>&1`
    if $?.to_i == 0
    # $compiler = "ftn95"
      $linker = "slink"
      $option = "/CFPP /DEFINE FTN95 1 /SILENT "
      $linker_option = ""
      $option_obj = "/binary"
      $ext_obj = "obj"
      $option_exe = "-out:"
      $dosish_path = true
      $gcov = false
		  $prof_genx = false
    end
  end
  #---------------------------------------------------

  $mpiexec_exist = false
  result_mpi = `#{where_or_which} mpiexec 2>&1`
  if $?.to_i == 0
    $mpiexec = "mpiexec"
    $mpiexec_exist = true
  end


  $linker = $compiler if !$linker
  $linker_option = $option if !$linker_option
  $option_obj = " -c -o " if !$option_obj
  $ext_obj    = "o"       if !$ext_obj
  $option_exe = " -o "    if !$option_exe

  result_which = `which ar 2>&1`
  $ar_ok = false
  $ar_ok = true if $?.to_i == 0
  $ar_ok = false if $compiler == "ftn95"
end
