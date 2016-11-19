
# Copyright (c) 2005-2010, 2012-2013, Andrew Hang Chen and contributors,
# All rights reserved.
# Licensed under the 3-clause BSD license.

require 'misc'

module RakeBaseDeps
  require 'rubygems'
  require 'fruit_processor'
  require 'rake/clean'
  require 'fileutils'

  include Rake::DSL if defined?(Rake::DSL)

  anchor_file_name='ROOT_ANCHOR'

  extensions = ["f90", "f95", "f03", "f08"]

  def opt_of_build_dir(build_dir)
    flag = build_dir
    if flag.size == 0 or flag == "./"
      flag = ""
    else
      flag = '"-I' + flag + '"' 
    end
    return flag
  end
  module_function :opt_of_build_dir

  def conv_dosish(target, source, if_dosish_path)
    target2 = ""
    source2 = ""
    if (if_dosish_path)
      target2 = target.gsub(%r{/}) { "\\" }
      source2 = source.gsub(%r{/}) { "\\" }
    else
      target2 = target
      source2 = source
    end
    return target2, source2
  end
  module_function :conv_dosish

  def copy_mod_to_build_dir(mods)
    mods.each do |module_file|
      if $build_dir != "" and $build_dir !~ /^\.\/*$/
        os_install File.expand_path(module_file), $build_dir
      end
    end
  end
  module_function :copy_mod_to_build_dir

  def coverage?(basename)
    if $coverage_fruit_f90 
      if basename == "fruit\." or basename == "fruit_util\."
        return true
      end
    end
    return false if basename =~ /_test\.$/
    return false if basename == "fruit\."
    return false if basename == "fruit_util\."
    return false if basename == "fruit_basket_gen\."
    return false if basename == "fruit_driver_gen\."
    return false if basename == "fruit_driver_mpi_gen\."
    return true
  end
  module_function :coverage?
#-------------------------------------------------------


  $goal = '' if !$goal

  $base_dir  = FruitProcessor.new.base_dir  if ! $base_dir
  $build_dir = FruitProcessor.new.build_dir if ! $build_dir

  $lib_bases  = {} if !$lib_bases
  $inc_dirs   = [] if !$inc_dirs
  $source_dir = "" if !$source_dir

  if not defined?(OBJ)
    # if rake_estimate.rb was not used, generate SRC and OBJ here
    SRC = FileList[]
    extensions.each{|fxx|
      SRC.concat(FileList['*.' + fxx])
    }
    SRC.sort!
    OBJ = SRC.ext($ext_obj)

    # generated files must be built last
    if OBJ.include?('fruit_basket_gen.' + $ext_obj)
      file 'fruit_basket_gen.' + $ext_obj =>  OBJ - ['fruit_basket_gen.' + $ext_obj,
                                                     'fruit_driver_gen.' + $ext_obj, 
                                                     'fruit_driver_mpi_gen.' + $ext_obj, ]
      file 'fruit_driver_gen.'     + $ext_obj =>  'fruit_basket_gen.' + $ext_obj
      file 'fruit_driver_mpi_gen.' + $ext_obj =>  'fruit_basket_gen.' + $ext_obj
    end
  end

  puts "rake_base_deps.rb: OBJ = " + OBJ.to_s if $show_info

  #-------------v
  # A_test.o depends on A.o if A.fxx exists
  OBJ.each{|a_obj|
    a = a_obj.to_s
    b = a.sub(/_test\.#{$ext_obj}$/, "\.#{$ext_obj}")
    if a != b then
      extensions.each{|fxx|
        b_src = b.sub(/\.#{$ext_obj}$/, "\.#{fxx}")
        if File.exist?(b_src)
          puts "rake_base_deps.rb: \"" + a + "\" depends on \"" + b + "\"" if $show_info
          file a => b
        end
      }
    end
  }
  #-------------^

  #------ coverage ------
  $for_coverage = []
  OBJ.each{|a_obj|
    basename = File.basename(a_obj, ".*") + '.'
    if coverage?(basename)
      #for gfortran + gcov
      if $gcov
        $for_coverage.push(basename + "gcda")
      end
      
      #for ifort
      if $prof_genx
        $for_coverage.push(basename + "dyn")
      end
    end
  }
  puts "For coverage: " + $for_coverage.to_s if $show_info

  rule ".gcda" => "." + $ext_obj
  rule ".dyn"  => "." + $ext_obj
  #------ coverage ------

  CLEAN.include([
    '*.o', '*.obj', '*.a', '*.mod',
    '*.gcda', '*.gcno', '*.gcov',
    '*.SPI', '*.SPL', '*.dpi', '*.dyn', 'CODE_COVERAGE.HTML',
    'fruit_*_gen.f90', '*fruit_driver', 'result*.xml',
    FruitProcessor.new.module_files(SRC, $build_dir)
  ])
  CLOBBER.include("#{$build_dir}/#{$goal}")

  task :default => [:deploy]

  # final goal link is depending on the libraries
  # file $goal => FruitProcessor.new.lib_base_files($lib_bases)
  # This is to resolve the .so files in LD_LIBRARY_PATH,
  # so if the .a file is not in the path, assume it is in one of the LD_LIBRARY_PATH
  FruitProcessor.new.lib_base_files($lib_bases).each do |lib|
    file $goal => lib if File.exist? lib
  end


  if !$source_dirs
    extensions.each{|fxx|
      rule '.' + $ext_obj => $source_dir + '%X.' + fxx do |t|
        Rake::Task[:dirs].invoke if Rake::Task.task_defined?('dirs')

        (name, source) = conv_dosish(t.name, t.source, $dosish_path)
        flag = opt_of_build_dir($build_dir)

        basename = File.basename(t.source, ".*") + '.'
        if $gcov and coverage?(basename)
          flag = flag + " " + $gcov
        end
        if $prof_genx and coverage?(basename)
          flag = flag + " " + $prof_genx
        end

        sh "#{$compiler} #{$option} #{$option_obj} " + 
           "#{name} #{source} #{flag} " + 
           "#{FruitProcessor.new.inc_flag($inc_dirs)}"

        copy_mod_to_build_dir(FileList["*.mod"])
      end
    }
  end

  if $source_dirs
    puts "======\ndependencies given by " + Pathname(__FILE__).basename.to_s + ":" if $show_info 
    puts "directries " + $source_dirs.to_s if $show_info
    $source_dirs.each{|dir|

      dir2 = Pathname.new(dir).cleanpath.to_s
      dir2 += "/" if dir2 != ""

      puts "in directry " + dir2 if $show_info

      extensions.each{|fxx|
        FileList[ dir2 + "*." + fxx ].each do |ff|
          source_fxx = Pathname.new(ff).cleanpath.to_s
          basename   = Pathname.new(ff).basename(fxx).to_s
          basename_o = basename + $ext_obj

          puts basename_o.to_s + " => " + source_fxx if $show_info 

          file basename_o => source_fxx do |t|
            Rake::Task[:dirs].invoke if Rake::Task.task_defined?('dirs')

            flag = opt_of_build_dir($build_dir)
            (name, source) = conv_dosish(t.name, source_fxx, $dosish_path)

            if $gcov and coverage?(basename)
              flag = flag + " " + $gcov
            end
            if $prof_genx and coverage?(basename)
              flag = flag + " " + $prof_genx
            end

            sh "#{$compiler} #{$option} #{$option_obj} " + 
               "#{name} #{source} #{flag} " + 
               "#{FruitProcessor.new.inc_flag($inc_dirs)}"

            copy_mod_to_build_dir(FileList["*.mod"])
          end
        end
      }
    }
    puts "======" if $show_info 
  end

  #`where ...` works on windows vista, 7 and 8. Not works on Windows XP.
  test_where = `where where 2>&1`
  if $?.to_i == 0
    where_or_which = "where"
  else
    where_or_which = "which"
  end

  file $goal => OBJ do
    if OBJ.size == 0
    elsif $goal =~ /.a$/
      result = `#{where_or_which} ar 2>&1`
      if $?.to_i == 0
        sh "ar cr #{$goal} #{OBJ}"
      end
    else
      lib_name_flag = FruitProcessor.new.lib_name_flag($lib_bases, $build_dir)

      lib_dir = FruitProcessor.new.lib_dir_flag($lib_bases, $build_dir)
      lib_dir = '"' + lib_dir + '"' if lib_dir.size > 0

      flag = $build_dir
      if flag.size == 0 or flag == "./"
        flag = ""
      else
        flag = '"-I' + flag + '"' 
      end

if $gcov
  flag = flag + " " + $gcov
end

      sh "#{$linker} #{$linker_option} #{flag} #{$option_exe}#{$goal} #{OBJ} #{lib_name_flag} #{lib_dir}"
    end
  end

  # generate directories
  task :dirs do
    if $build_dir != "" and $build_dir !~ /^\.\/*$/
      Dir.mkdir $build_dir unless File.exist?($build_dir)
    end
  end

  task :deploy => $goal do
    if $goal.length > 0
      if $build_dir != "" and $build_dir !~ /^\.\/*$/
        os_install "#{Dir.pwd}/#{$goal}", "#{$build_dir}/#{$goal}"
      end
    end
  end

  task :anchor_root do
    FileUtils.touch anchor_file_name if not File.exist? anchor_file_name
  end

  task :gen do
    FruitProcessor.new.pre_process
  end

  task :spec do
    FruitProcessor.new.spec_report
  end

end
