%w(colorize fileutils rake/clean).each do |gem|
  begin
    require gem
  rescue LoadError
    warn "Install the #{gem} gem:\n $ (sudo) gem install #{gem}".magenta
    exit 1
  end
end
require 'shellwords'

CLEAN.clear_exclude.exclude { |fn| fn.pathmap("%f").downcase == "core" }

# Configurazione esterna (opzionale)
if File.exist?(File.expand_path('../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../Rakefile_configure.rb'
elsif File.exist?(File.expand_path('../../Rakefile_configure.rb', File.dirname(__FILE__)))
  require_relative '../../Rakefile_configure.rb'
else
  COMPILE_DEBUG      = false
  COMPILE_DYNAMIC    = false
  COMPILE_EXECUTABLE = true
end

# ----------------------------------------------------------------------------
# Rilevamento OS
# ----------------------------------------------------------------------------
case RUBY_PLATFORM
when /darwin/
  OS = :mac
when /linux|cygwin/
  OS = :linux
when /msys/
  OS = :mingw
else
  OS = :win
end

def build_type
  COMPILE_DEBUG ? 'Debug' : 'Release'
end

def install_prefix
  File.expand_path('lib', __dir__)
end

def project_root
  File.expand_path(__dir__)
end

def build_dir
  File.join(project_root, 'build')
end

def cmake_generator
  'Ninja'
end

def platform_suffix
  case OS
  when :mac
    'osx'
  when :linux
    'linux'
  when :mingw
    'mingw'
  else
    'win'
  end
end

def cleanup_install_layout
  root_lib = File.join(project_root, 'lib')
  public_lib = File.join(root_lib, 'lib')
  FileUtils.rm_rf File.join(root_lib, 'Users')

  # Root-level libraries are invalid for this project. The real Astro archive
  # and its generic symlink must live in lib/lib only.
  Dir.glob(File.join(root_lib, 'lib*.{a,dylib,so}'), File::FNM_EXTGLOB).each do |path|
    FileUtils.rm_rf(path) if File.file?(path) || File.symlink?(path)
  end

  # Dependency artifacts belong to lib3rd/lib, not to lib/lib. Keep only Astro
  # artifacts in the public library directory.
  dependency_patterns = %w[
    libGenericContainer* libLua* libquarticRootsFlocke* libUtilsLite*
  ]
  dependency_patterns.each do |pattern|
    Dir.glob(File.join(public_lib, pattern)).each do |path|
      FileUtils.rm_rf(path) if File.file?(path) || File.symlink?(path)
    end
  end
end

def cmake_configure_command(build_tests: false)
  args = [
    "-G", cmake_generator,
    "-B", build_dir,
    "-DCMAKE_BUILD_TYPE=#{build_type}",
    "-DCMAKE_INSTALL_PREFIX=#{install_prefix}",
    "-DBUILD_SHARED_LIBS=#{COMPILE_DYNAMIC ? 'ON' : 'OFF'}",
    "-DBUILD_TESTING=#{build_tests ? 'ON' : 'OFF'}",
    "-DASTRO_INSTALL=ON",
    "-DASTRO_BUILD_BENCHMARKS=OFF",
    "-DASTRO_STRICT_WARNINGS=OFF",
    project_root
  ]
  args.join(' ')
end

def cmake_build_command(target = nil)
  cmd = ["--build", build_dir, "--config", build_type]
  cmd += ["--target", target] if target
  cmd += ["--parallel"]
  cmd.join(' ')
end

# ----------------------------------------------------------------------------
# Task principali
# ----------------------------------------------------------------------------
desc "default task --> build"
task :default => :build

desc "git clean reset"
task :git_clean do
  sh "git reset --hard"
  sh "git clean -d -x -f"
end

desc "Configure CMake (without building)"
task :configure do
  puts "Configuring CMake in #{build_dir}...".green
  sh "cmake " + cmake_configure_command(build_tests: false)
end

desc "Build and install the library (does NOT compile tests)"
task :build => :configure do
  puts "Compiling and installing library...".green
  sh "cmake " + cmake_build_command('install')
  cleanup_install_layout
end

desc "Configure CMake with tests enabled"
task :configure_tests do
  puts "Configuring CMake with tests in #{build_dir}...".green
  sh "cmake " + cmake_configure_command(build_tests: true)
end

desc "Compile test executables"
task :build_tests => :configure_tests do
  puts "Compiling tests...".green
  sh "cmake " + cmake_build_command('Astro_all_tests')
  cleanup_install_layout
end

desc "Run all tests (compiles tests if needed)"
task :test => :build_tests do
  puts "Running tests...".green
  sh "ctest --test-dir #{build_dir.shellescape} --build-config #{build_type.shellescape} --output-on-failure"
end

desc "Build and run all tests (alias for test)"
task :run => :test do
  # same as test
end

desc "Clean build artifacts (keeps lib/ and lib3rd/)"
task :clean_build do
  puts "Removing build directory only...".green
  FileUtils.rm_rf build_dir
end

desc "Full clean (removes build, lib, lib3rd, bin, object files)"
task :clean do
  puts "Full cleaning...".green
  FileUtils.rm_rf build_dir
  FileUtils.rm_rf File.join(project_root, 'lib')
  FileUtils.rm_rf File.join(project_root, 'lib3rd')
  FileUtils.rm_rf File.join(project_root, 'bin')
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.o'))
  FileUtils.rm_f Dir.glob(File.join(project_root, '**', '*.obj'))
end

desc "Package using CPack"
task :cpack do
  puts "Creating packages...".green
  Dir.chdir(build_dir) do
    sh "cpack -C #{build_type} CPackConfig.cmake"
    sh "cpack -C #{build_type} CPackSourceConfig.cmake"
  end
end

CLOBBER.include []
