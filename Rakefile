require_relative "./cmake_utils/Rakefile_common.rb"

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win, [:year, :bits] do |t, args|

  args.with_defaults( :year => "2017", :bits => "x64" )

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd_cmake = cmake_generation_command(args.bits,args.year) + cmd_cmake_build()

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

desc "compile for OSX/LINUX/MINGW"
task :build_common do

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  cmd_cmake = "cmake " + cmd_cmake_build()

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
task :build_linux => :build_common do
end

desc "compile for OSX"
task :build_osx => :build_common do
end

desc "compile for MINGW"
task :build_mingw => :build_common do
end
