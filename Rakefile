%w(colorize rake fileutils).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end

require 'rake/clean'

CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux/
  OS = :linux
when /cygwin|mswin|mingw|bccwin|wince|emx/
  OS = :win
when /msys/
  OS = :win
end

require_relative "./Rakefile_common.rb"

file_base = File.expand_path(File.dirname(__FILE__)).to_s

cmd_cmake_build = ""
if COMPILE_EXECUTABLE then
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=ON '
else
  cmd_cmake_build += ' -DEB_ENABLE_TESTS:VAR=OFF '
end
if COMPILE_DYNAMIC then
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=ON '
else
  cmd_cmake_build += ' -DEB_BUILD_SHARED:VAR=OFF '
end
if COMPILE_DEBUG then
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Debug --loglevel=STATUS '
else
  cmd_cmake_build += ' -DCMAKE_BUILD_TYPE:VAR=Release --loglevel=STATUS '
end

desc "default task --> build"
task :default => :build

FileUtils.cp './cmake/CMakeLists-cflags.txt', './submodules/Utils/cmake/CMakeLists-cflags.txt'
FileUtils.cp './cmake/CMakeLists-cflags.txt', './submodules/quarticRootsFlocke/cmake/CMakeLists-cflags.txt'
FileUtils.cp './cmake/CMakeLists-cflags.txt', './submodules/GenericContainer/CMakeLists-cflags.txt'

desc "run tests"
task :run do
  puts "ASTRO run tests".green
  case OS
  when :mac,:linux
    Dir.glob('./bin/test_*').each do |exe|
      next unless File.exist?(exe)
      puts "execute #{exe}".yellow
      sh exe
    end
  when :win
    Dir.glob('./bin/test_*.exe').each do |exe|
      next unless File.exist?(exe)
      puts "execute #{exe}".yellow
      system(exe)
    end
  end
end

desc "build ASTRO"
task :build do
  case OS
  when :mac
    puts "ASTRO build (osx)".green
    Rake::Task[:build_osx].invoke
  when :linux
    puts "ASTRO build (linux)".green
    Rake::Task[:build_linux].invoke
  when :win
    puts "ASTRO build (windows)".green
    Rake::Task[:build_win].invoke
  end
end


desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  dir = "vs_#{args.year}_#{args.bits}"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = win_vs(args.bits,args.year) + cmd_cmake_build

  puts "run CMAKE for ASTRO".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for ASTRO".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for OSX"
task :build_common, [:os] do |t, args|

  args.with_defaults( :os => "osx" )

  dir = "build"

  FileUtils.rm_rf   dir
  FileUtils.mkdir_p dir
  FileUtils.cd      dir

  cmd_cmake = "cmake " + cmd_cmake_build

  puts "run CMAKE for ASTRO".yellow
  sh cmd_cmake + ' ..'
  puts "compile with CMAKE for ASTRO".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL+QUIET
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL+QUIET
  end

  FileUtils.cd '..'
end

desc "compile for LINUX"
task :build_linux do
  Rake::Task[:build_common].invoke("linux")
end

desc "compile for OSX"
task :build_osx do
  Rake::Task[:build_common].invoke("osx")
end
