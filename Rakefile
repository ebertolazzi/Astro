require_relative "./cmake_utils/Rakefile_common.rb"

desc "compile for Visual Studio [default year=2017, bits=x64]"
task :build_win do
  # check architecture
  case `where cl.exe`.chop
  when /x64\\cl\.exe/
    VS_ARCH = 'x64'
  when /amd64\\cl\.exe/
    VS_ARCH = 'x64'
  when /bin\\cl\.exe/
    VS_ARCH = 'x86'
  else
    raise RuntimeError, "Cannot determine architecture for Visual Studio".red
  end

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "run CMAKE for ASTRO".yellow
  sh "cmake -G Ninja -DBITS:VAR=#{VS_ARCH} " + cmd_cmake_build() + ' ..'

  puts "compile with CMAKE for ASTRO".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake --build . --config Release --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

desc "compile for OSX/LINUX/MINGW"
task :build_common do

  FileUtils.rm_rf   "build"
  FileUtils.mkdir_p "build"
  FileUtils.cd      "build"

  puts "run CMAKE for ASTRO".yellow
  sh "cmake -G Ninja " + cmd_cmake_build() + ' ..'

  puts "compile with CMAKE for ASTRO".yellow
  if COMPILE_DEBUG then
    sh 'cmake --build . --config Debug --target install '+PARALLEL
  else
    sh 'cmake  --build . --config Release  --target install '+PARALLEL
  end

  FileUtils.cd '..'
end

task :build_linux => :build_common do end
task :build_osx   => :build_common do end
task :build_mingw => :build_common do end
